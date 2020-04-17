-- Gkeyll
-- Two-fluid 5m simulation
-- Current filamentation instability (Weibel)
-- Petr

----------------------------------------------------------------------
-- Simulation parameters ---------------------------------------------
pi = Lucee.Pi
-- Constatnts
gasGamma = 5./3.
chargeElc = -1.0
massElc = 1.0
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1.0/math.sqrt(epsilon0*mu0)

-- Initial conditions
nElc10 = 0.5
nElc20 = 0.5
uxElc10 = 0.0
uyElc10 = 0.3
uzElc10 = 0.0
uxElc20 = 0.0
uyElc20 = -0.3
uzElc20 = 0.0
TElc10 = 0.01
TElc20 = 0.01
k0 = 0.4
perturb = 1e-3
-- IC automatically calculated
vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)

-- Domain and time
numCells = {384} -- Nx
lowerBoundary = {0.0} -- xLow
upperBoundary = {2*pi/k0} -- xUp
periodicDirs = {0}
cfl = 0.9
tStart = 0.0
tEnd = 150.0
numFrames = 300

-- Other setup
numFluids = 2
Lucee.IsRestarting = false
Lucee.RestartFrame = -1
limiter = "van-leer"
elcErrorSpeedFactor = 0
mgnErrorSpeedFactor = 1


----------------------------------------------------------------------
-- Convenience and Debugging -----------------------------------------
log = function(...) Lucee.logInfo(string.format(...)) end

function HERE()
   local info = debug.getinfo(2)
   str = string.format("HERE: %d %s: %s", info.currentline,
		       info.source, tostring(info.name))
   Lucee.logInfo(string.format(str))
end


----------------------------------------------------------------------
-- Setup verification ------------------------------------------------
log("Setup verification:")
log(" * tEnd = %g,  numFrames = %d", tEnd, numFrames)
log(" * Speed of light = %g", lightSpeed)
log(" * Vthe1/c = %g", vthElc10/lightSpeed)
log(" * Vthe2/c = %g", vthElc20/lightSpeed)


----------------------------------------------------------------------
-- Computational domain, Data fields ---------------------------------
decomp = DecompRegionCalc1D.CartGeneral {}
grid = Grid.RectCart1D {
   lower = lowerBoundary,
   upper = upperBoundary,
   cells = numCells,
   decomposition = decomp,
   periodicDirs = periodicDirs,
}

createField = function(numComponents)
   return DataStruct.Field1D {
      onGrid = grid,
      numComponents = numComponents,
      ghost = {2, 2},
   }
end
q = createField(18)
qNew = createField(18)
qDup = createField(18)

-- aliases to various sub-systems
getFields = function(fld)
   return fld:alias(0, 5), fld:alias(5, 10), fld:alias(10, 18)
end
elcFluid1, elcFluid2, emField = getFields(q)
elcFluid1New, elcFluid2New, emFieldNew = getFields(qNew)


----------------------------------------------------------------------
-- Boundary conditions ----------------------------------------------- 
function applyBc(fld, tCurr, tAdv)
   --fld:applyCopyBc(0, "lower")
   --fld:applyCopyBc(0, "upper")
   fld:sync()
end


----------------------------------------------------------------------
-- Initial conditions setup ------------------------------------------
function init(x, y, z)
   local rho1 = massElc*nElc10
   local mx1 = rho1*uxElc10
   local my1 = rho1*uyElc10
   local mz1 = rho1*uzElc10
   local p1 = nElc10*TElc10/(gasGamma-1) +
      0.5*(mx1*mx1 + my1*my1 + mz1*mz1)/rho1


   local rho2 = massElc*nElc20
   local mx2 = rho2*uxElc20
   local my2 = rho2*uyElc20
   local mz2 = rho2*uzElc20
   local p2 = nElc20*TElc20/(gasGamma-1) +
      0.5*(mx2*mx2 + my2*my2 + mz2*mz2)/rho2

   local Ex = 0.0
   local Ey = 0.0
   local Ez = 0.0
   local Bx = 0.0
   local By = 0.0
   local Bz = perturb*math.sin(k0*x)

   return rho1, mx1, my1, mz1, p1, rho2, mx2, my2, mz2, p2, Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
end


----------------------------------------------------------------------
-- Equation solvers --------------------------------------------------
eulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
eulerEqnLax = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
}

createSlvr = function(eqn, input, output, dir, limiter)
   local slvr = Updater.WavePropagation1D {
      onGrid = grid,
      equation = eqn,
      limiter = limiter,
      cfl = cfl,
      cflm = 1.1*cfl,
      updateDirections = {dir}
   }
   slvr:setIn( {input} )
   slvr:setOut( {output} )
   return slvr
end

elc1SlvrDir0 = createSlvr(eulerEqn, elcFluid1, elcFluid1New, 0, limiter)
elc2SlvrDir0 = createSlvr(eulerEqn, elcFluid2, elcFluid2New, 0, limiter)
emSlvrDir0 = createSlvr(maxwellEqn, emField, emFieldNew, 0, limiter)

elc1SlvrDir0Lax = createSlvr(eulerEqnLax, elcFluid1, elcFluid1New, 0, "zero")
elc2SlvrDir0Lax = createSlvr(eulerEqnLax, elcFluid2, elcFluid2New, 0, "zero")
emSlvrDir0Lax = createSlvr(maxwellEqn, emField, emFieldNew, 0, "zero")

slvrs = {
   {elc1SlvrDir0, elc2SlvrDir0, emSlvrDir0},
}
slvrsLax = {
   {elc1SlvrDir0Lax, elc2SlvrDir0Lax, emSlvrDir0Lax},
}


----------------------------------------------------------------------
-- Source solvers ----------------------------------------------------
sourceSlvr = Updater.ImplicitFiveMomentSrc1D {
   onGrid = grid,
   numFluids = 2,
   charge = {chargeElc, chargeElc},
   mass = {massElc, massElc},
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

function updateSource(elc1In, elc2In, emIn, tCurr, tAdv)
   sourceSlvr:setOut({elc1In, elc2In, emIn})
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(tAdv)
end


----------------------------------------------------------------------
-- Field updaters ----------------------------------------------------
function updateFluidsAndField(tCurr, tAdv)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tAdv - tCurr)
   local useLaxSolver = False

   for i, slvr in ipairs(slvrs[1]) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(tAdv)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((eulerEqn:checkInvariantDomain(elcFluid1New) == false)
    or (eulerEqn:checkInvariantDomain(elcFluid2New) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   return myStatus, myDtSuggested, useLaxSolver
end

function updateFluidsAndFieldLax(tCurr, tAdv)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tAdv - tCurr)
   for i, slvr in ipairs(slvrsLax[1]) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(tAdv)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

function updateSystem(tCurr, tAdv)
   local dtHalf = 0.5*(tAdv - tCurr)

   updateSource(elcFluid1, elcFluid2, emField, tCurr, tCurr + dtHalf)
   applyBc(q, tCurr, tAdv)

   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, tAdv)

   updateSource(elcFluid1New, elcFluid2New, emFieldNew, tCurr, tCurr + dtHalf)
   applyBc(qNew, tCurr, tAdv)

   return status, dtSuggested, useLaxSolver
end

function updateSystemLax(tCurr, tAdv)
   local dtHalf = 0.5*(tAdv - tCurr)

   updateSource(elcFluid1, elcFluid2, emField, tCurr, tCurr + dtHalf)
   applyBc(q, tCurr, tAdv)

   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, tAdv)

   updateSource(elcFluid1New, elcFluid2New, emFieldNew, tCurr, tCurr + dtHalf)
   applyBc(qNew, tCurr, tAdv)

   return status, dtSuggested
end


----------------------------------------------------------------------
-- Diagnostics and Output --------------------------------------------
emEnergy = DataStruct.DynVector {numComponents = 1}
emEnergyCalc = Updater.IntegrateField1D {
   onGrid = grid,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*(ex*ex + ey*ey + ez*ez) +
	 0.5/mu0*(bx*bx + by*by + bz*bz)
   end,
}
emEnergyCalc:setIn({emField})
emEnergyCalc:setOut({emEnergy})

ExEnergy = DataStruct.DynVector {numComponents = 1}
ExEnergyCalc = Updater.IntegrateField1D {
   onGrid = grid,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*ex*ex
	       end,
}
ExEnergyCalc:setIn({emField})
ExEnergyCalc:setOut({ExEnergy})

EyEnergy = DataStruct.DynVector {numComponents = 1}
EyEnergyCalc = Updater.IntegrateField1D {
   onGrid = grid,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5*epsilon0*ey*ey
	       end,
}
EyEnergyCalc:setIn({emField})
EyEnergyCalc:setOut({EyEnergy})

BzEnergy = DataStruct.DynVector {numComponents = 1}
BzEnergyCalc = Updater.IntegrateField1D {
   onGrid = grid,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return 0.5/mu0*bz*bz
	       end,
}
BzEnergyCalc:setIn({emField})
BzEnergyCalc:setOut({BzEnergy})

totalPtclEnergy = DataStruct.DynVector {numComponents = 1}
totalPtclEnergyCalc = Updater.IntegrateField1D {
   onGrid = grid,
   integrand = function (rho1, mx1, my1, mz1, p1, rho2, mx2, my2, mz2, p2, ex, ey, ez, bx, by, bz, e1, e2)
		  return p1 + p2
	       end,
}
totalPtclEnergyCalc:setIn({qNew} )
totalPtclEnergyCalc:setOut({totalPtclEnergy})

function calcDiagnostics(tCurr, tAdv)
   for i, diag in ipairs({emEnergyCalc, ExEnergyCalc, EyEnergyCalc,
			  BzEnergyCalc, totalPtclEnergyCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tAdv)
   end
end

function writeFields(frame, tCurr)
   qNew:write( string.format("q_%03d.h5", frame), tCurr)

   emEnergy:write(string.format("emEnergy_%03d.h5", frame))
   ExEnergy:write(string.format("ExEnergy_%03d.h5", frame))
   EyEnergy:write(string.format("EyEnergy_%03d.h5", frame))
   BzEnergy:write(string.format("BzEnergy_%03d.h5", frame))
   totalPtclEnergy:write(string.format("totalPtclEnergy_%03d.h5", frame))
end

function writeDead(tCurr)
   q:write( string.format("q_dead.h5"), tCurr)
   qNew:write( string.format("qNew_dead.h5"), tCurr)
end


----------------------------------------------------------------------
-- Main function -----------------------------------------------------
function runSimulation(tStart, tEnd, numFrames, initDt)
   local frame = 1
   local tFrame = (tEnd - tStart)/numFrames
   local tNextFrame = tFrame
   local step = 1
   local stepFrame = 1
   local tCurr = tStart
   local dt = initDt
   local status = true
   local dtSuggested = initDt
   local useLaxSolver = false
   
   if Lucee.IsRestarting then
      fileName = "q_" .. Lucee.RestartFrame .. ".h5"
      tCurr = q:read(fileName)
      if not tCurr then
         tCurr = tStart + (tEnd - tStart) * Lucee.RestartFrame / numFrames
      end
      frame = Lucee.RestartFrame + 1
      tNextFrame = tCurr + tFrame
      log('\nRestarting from frame %d tCurr = %g\n', Lucee.RestartFrame, tCurr)
   else
      q:set(init)
   end
   applyBc(q)
   qNew:copy(q)
   qNew:sync()
   
   if not Lucee.IsRestarting then
      calcDiagnostics(tStart, tStart)
      writeFields(0, tStart)
   end
   
   while true do
      qDup:copy(q)
      
      if (tCurr + dt > tEnd) then
	 dt = tEnd - tCurr
      end
      
      -- advance fluids and fields
      if useLaxFlux then
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g; using Lax fluxes", step, stepFrame, frame, tCurr, dt)
         status, dtSuggested = updateSystemLax(tCurr, tCurr + dt)
         if status then
            useLaxFlux = false
         end
      else
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g",
	     step, stepFrame, frame, tCurr, dt)
         status, dtSuggested, useLaxFlux = updateSystem(tCurr, tCurr + dt)
      end
      
      if (status == false) then
	 -- time-step too large
	 log(" ** dt %g too large! Will retake with dt %g", dt, dtSuggested)
	 dt = dtSuggested
	 q:copy(qDup)
      elseif (useLaxSolver == true) then
	 -- negative density/pressure occured
	 log (" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr + dt)
	 q:copy(qDup)
      else
	 -- check if a nan occured
	 if (qNew:hasNan()) then
	    log(" ** NaN occured at %g! Stopping simulation", tCurr + dt)
	    writeDead(tCurr + dt)
	    break
	 end
	 
	 calcDiagnostics(tCurr, tCurr + dt)
	 q:copy(qNew)        
	 tCurr = tCurr + dt
	 dt = dtSuggested
	 step = step + 1
	 stepFrame = stepFrame + 1
	 
	 if (tCurr > tNextFrame or tCurr >= tEnd) then
	    log("Writing data at time %g (frame %d) ...\n", tCurr, frame)
	    writeFields(frame, tCurr)
	    frame = frame + 1
	    tNextFrame = tNextFrame + tFrame
	    stepFrame = 1
	 end  
	 
	 if (tCurr >= tEnd) then
           break
	 end
      end 
   end 
end

----------------------------------------------------------------------
-- The Thing! --------------------------------------------------------
runSimulation(tStart, tEnd, numFrames, tEnd - tStart)
