-- Vlasov-Maxwell solver: Units are SI

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- Permeabilty of free space
lightSpeed = 1.0 -- speed of light

Te_Ti = 1.0 -- ratio of electron to ion temperature
elcTemp = 0.01 -- electron temperature
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

ionTemp = elcTemp/Te_Ti -- ion temperature
ionMass = 1836 -- ion mass
ionCharge = 1.0 -- ion charge

-- thermal speeds
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)

-- simulation parameters
len = 1.0
delta = 0.25
uinf = 0.3
k0 = 1.4
n0 = 1.0
perturb = 1e-3

-- domain size and simulation time
xlower = 0.0 -- lower bounds of domain
xupper = 2*Lucee.Pi/k0 -- upper bounds of domain

LX = xupper-xlower

-- Resolution, time-stepping etc.
NX = 128
NVX = 32
NVY = 32
--NVZ = 16
tStart = 0.0 -- start time 
tEnd = 60
nFrames = 60

cfl = 0.5*(1.0/3.0)/(2*polyOrder+1)

-- compute max thermal speed to set velocity space extents
VXL_ELC, VXU_ELC = -0.80, 0.80 -- -10.0*vtElc, 10.0*vtElc
VYL_ELC, VYU_ELC = -0.80, 0.80 -- -10.0*vtElc, 10.0*vtElc
VL_ION, VU_ION = -6.0*vtIon, 6.0*vtIon

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames=%d", tEnd, nFrames))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Cell size=%g", LX/NX))
log(string.format("Light speed=%g", lightSpeed))
log(string.format("Time-step from light speed=%g", cfl*LX/NX/lightSpeed))
log(string.format("Electron thermal speed=%g", vtElc))
log(string.format("Ion thermal speed=%g", vtIon))
log(string.format("Configuration domain extents = [%g,%g]", -LX/2, LX/2))
--log(string.format("Electron velocity domain extents = [%g,%g]", VL_ELC, VU_ELC))
log(string.format("Ion velociy domain extents = [%g,%g]", VL_ION, VU_ION))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc3D.CartProd { cuts = {2,2,2} }
confDecomp = DecompRegionCalc1D.SubCartProd3D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid for electrons
phaseGridElc = Grid.RectCart3D {
   lower = {xlower, VXL_ELC, VYL_ELC},
   upper = {xupper, VXU_ELC, VYU_ELC},
   cells = {NX, NVX, NVY},
   decomposition = phaseDecomp,
   periodicDirs = {0},   
}
-- phase space grid for ions
phaseGridIon = Grid.RectCart3D {
   lower = {xlower, VL_ION, VL_ION},
   upper = {xupper, VU_ION, VU_ION},
   cells = {NX, NVX, NVY},
   decomposition = phaseDecomp,
   periodicDirs = {0},
}

-- configuration space grid (same for electrons and ions)
confGrid = Grid.RectCart1D {
   lower = {xlower},
   upper = {xupper},
   cells = {NX},
   decomposition = confDecomp,
   periodicDirs = {0},
}
-- phase-space basis functions for electrons and ions
phaseBasisElc = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridElc,
   polyOrder = polyOrder,
}
phaseBasisIon = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridIon,
   polyOrder = polyOrder,
}
-- configuration-space basis functions (shared by both species)
confBasis = NodalFiniteElement1D.LagrangeTensor {
   onGrid = confGrid,
   polyOrder = polyOrder,
   nodeLocation = "uniform",
}

-- distribution function for electrons
distfElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- extra fields for performing RK update
distfNewElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- Electron number density
numDensityElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion number density
numDensityIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron momentum
momentumElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}

-- net current
current = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for adding to EM fields (this perhaps is not the best way to do it)
emSource = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}

-- EM field
em = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
em1 = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}
emNew = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 8*confBasis:numNodes(),
   ghost = {1, 1},
}

--------------------------------
-- INITIAL CONDITION UPDATERS --
--------------------------------

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(vdrift_x, vdrift_y, vt, vx, vy)
   local v2 = (vx-vdrift_x)^2 + (vy-vdrift_y)^2
   --local norm = math.sqrt(2*Lucee.Pi*vt^2)
   return n0/(2*Lucee.Pi*vt^2)*math.exp(-v2/(2*vt^2))
end

-- get the initial perturbation
function perturbation(x)
   return perturb*math.sin(k0*x)
end

-- get the initial drift velocity 
function drift(x)
   return uinf+0.5*delta*math.tanh(x/len)
end

-- updater to initialize electron distribution function
initDistfElc = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,vx,vy,t)
      --local vD = drift(x)
      return 0.5*maxwellian( 0.0, uinf, vtElc, vx, vy) +
	 0.5*maxwellian(0.0, -uinf, vtElc, vx, vy) --+ 
	 --perturbation(x)/(12*vtElc)^2
   end
}

-- updater to initialize ion distribution function
initDistfIon = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,vx,vy,t)
      return maxwellian(0.0, 0.0, vtIon, vx, vy)
   end
}

-- updater to initialize EM fields
initField = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      -- no fields initially
      return 0.0, 0.0, 0.0, 0.0, 0.0, perturbation(x), 0.0, 0.0
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.NodalVlasov1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   cfl = cfl,
   charge = elcCharge,
   mass = elcMass,
}
vlasovSolverIon = Updater.NodalVlasov1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   cfl = cfl,
   charge = ionCharge,
   mass = ionMass,
}

-- Maxwell equation object
maxwellEqn = HyperEquation.PhMaxwell {
   -- speed of light
   lightSpeed = lightSpeed,
   -- factor for electric field correction potential speed
   elcErrorSpeedFactor = 0.0,
   -- factor for magnetic field correction potential speed
   mgnErrorSpeedFactor = 1.0,
   -- numerical flux to use: one of "upwind" or "central"
   numericalFlux = "upwind",
}

-- updater to solve Maxwell equations
maxwellSlvr = Updater.NodalDgHyper1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- equation system to solver
   equation = maxwellEqn,
   -- CFL number
   cfl = cfl,
}


-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,   
   moment = 0,
}
numDensityCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   moment = 0,
}

-- Updater to compute momentum
momentumCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 1,
}

-- This strange looking updater copies the currents into the EM source
-- field. Perhaps this is not the best way to do things, and one can
-- imagine a source updater which adds current sources to the dE/dt
-- Maxwell equation
copyToEmSource = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {0, 1},
   targetComponents = {0, 1},
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields


-- apply boundary conditions to distribution functions
function applyDistFuncBc(curr, dt, fldElc, fldIon)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fldElc})
   end
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {fldIon})
   end

   for i,fld in ipairs({fldElc, fldIon}) do
   end

   -- sync the distribution function across processors
   fldElc:sync()
   fldIon:sync()
end

-- apply BCs to EM fields
function applyEmBc(curr, dt, EM)
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {EM})
   end
   for i,bc in ipairs({}) do
      runUpdater(bc, curr, dt, {}, {EM})
   end

   -- sync EM fields across processors
   EM:sync()
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

emEnergy = DataStruct.DynVector { numComponents = 1 }
emEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*(ex^2+ey^2+ez^2) + 0.5/mu0*(bx^2+by^2+bz^2)
   end,
}
emEnergyCalc:setIn( {em} )
emEnergyCalc:setOut( {emEnergy} )

Bz2 = DataStruct.DynVector { numComponents = 1 }
Bz2Calc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return bz^2
   end,
}
Bz2Calc:setIn( {em} )
Bz2Calc:setOut( {Bz2} )

----------------------
-- SOLVER UTILITIES --
----------------------

-- generic function to run an updater
function runUpdater(updater, currTime, timeStep, inpFlds, outFlds)
   updater:setCurrTime(currTime)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(currTime+timeStep)
end

-- function to calculate number density
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end
-- function to calculate momentum density
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end

-- function to compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
   -- number density
   calcNumDensity(numDensityCalcElc, curr, dt, distfElcIn, numDensityElc)
end

-- function to update Vlasov equation
function updateVlasovEqn(vlasovSlvr, curr, dt, distfIn, emIn, distfOut)
   return runUpdater(vlasovSlvr, curr, dt, {distfIn, emIn}, {distfOut})
end

-- solve maxwell equation
function updateMaxwellEqn(curr, dt, emIn, emOut)
   return runUpdater(maxwellSlvr, curr, dt, {emIn}, {emOut})
end

-- function to compute diagnostics
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({Bz2Calc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

----------------------------
-- Time-stepping routines --
----------------------------

-- take single RK step
function rkStage(tCurr, dt, elcIn, ionIn, emIn, elcOut, ionOut, emOut)
   -- update distribution functions and homogenous Maxwell equations
   local stElc, dtElc = updateVlasovEqn(vlasovSolverElc, tCurr, dt, 
					elcIn, emIn, elcOut)
   ionOut:copy(ionIn) -- not updating ion equations, 
                      -- so just copy old solution to output field
   local stEm, dtEm = updateMaxwellEqn(tCurr, dt, emIn, emOut)
   
   -- compute currents from old values
   calcMomentum(momentumCalcElc, tCurr, dt, elcIn, momentumElc)
   current:combine(elcCharge, momentumElc, ionCharge, momentumIon)
   
   -- copy into EM sources
   runUpdater(copyToEmSource, tCurr, dt, {current}, {emSource})
   -- add in current source to Maxwell equation output
   emOut:accumulate(-dt/epsilon0, emSource)
   
   if (stElc == false) or (stIon == false) or (stEm == false)  then
      return false, math.min(dtElc, dtEm)
   end
   return true, math.min(dtElc, dtEm)
end

function rk3(tCurr, myDt)
   local myStatus, myDtSuggested
   -- RK stage 1
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distfElc, distfIon, em, distf1Elc, distf1Ion, em1)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, dt, em1)

   -- RK stage 2
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Elc, distf1Ion, em1, distfNewElc, distfNewIon, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   em1:combine(3.0/4.0, em, 1.0/4.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, dt, em1)   

   -- RK stage 3
   myStatus, myDtSuggested = rkStage(tCurr, myDt, distf1Elc, distf1Ion, em1, distfNewElc, distfNewIon, emNew)
   if (myStatus == false)  then
      return false, myDtSuggested
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   em1:combine(1.0/3.0, em, 2.0/3.0, emNew)
   -- apply BC
   applyDistFuncBc(tCurr, dt, distf1Elc, distf1Ion)
   applyEmBc(tCurr, dt, em1)   

   distfElc:copy(distf1Elc)
   em:copy(em1)

   return true, myDtSuggested
end

-- make duplicates in case we need them
distfDupElc = distfElc:duplicate()
emDup = em:duplicate()

-- function to advance solution from tStart to tEnd
function runSimulation(tStart, tEnd, nFrames, initDt)
   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- the grand loop 
   while true do
      distfDupElc:copy(distfElc)
      emDup:copy(em)
      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      -- advance particles and fields
      log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)
      if (status == false) then
	 -- time-step too large
	 log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 distfElc:copy(distfDupElc)
	 em:copy(emDup)
      else
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)
	 -- write out data
	 if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then
	    log (string.format(" Writing data at time %g (frame %d) ...\n", tCurr+myDt, frame))
	    writeFields(frame, tCurr+myDt)
	    frame = frame + 1
	    nextIOt = nextIOt + tFrame
	    step = 0
	 end

	 tCurr = tCurr + myDt
	 myDt = dtSuggested
	 step = step + 1
	 -- check if done
	 if (tCurr >= tEnd) then
	    break
	 end
      end 
   end -- end of time-step loop
   return dtSuggested
end

-- Write out data frame 'frameNum' with at specified time 'tCurr'
function writeFields(frameNum, tCurr)
   -- distribution functions
   distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)

   -- EM field
   em:write(string.format("em_%d.h5", frameNum), tCurr)   

   -- compute moments and write them out
   calcMoments(tCurr, 0.0, distfElc, distfIon)
   numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
   momentumElc:write(string.format("momentumElc_%d.h5", frameNum), tCurr)
   ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)

   -- diagnostics
   Bz2:write( string.format("Bz2_%d.h5", frameNum) )
end


----------------------------
-- RUNNING THE SIMULATION --
----------------------------
-- apply initial conditions for electrons and ion
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
runUpdater(initField, 0.0, 0.0, {}, {em})

-- print IC timing information
log(string.format("Electron IC updater took = %g", 
		  initDistfElc:totalAdvanceTime()))
log(string.format("Ion IC updater took = %g", initDistfIon:totalAdvanceTime()))

-- apply BCs
applyDistFuncBc(0.0, 0.0, distfElc, distfIon)
applyEmBc(0.0, 0.0, em)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial fields
writeFields(0, 0.0)

-- ions are not evolved
calcNumDensity(numDensityCalcIon, 0, 0, distfIon, numDensityIon)
numDensityIon:write("numDensityIon.h5", 0)
calcMomentum(momentumCalcIon, 0, 0, distfIon, momentumIon)
momentumIon:write("momentumIon.h5", 0)

-- run the whole thing
initDt = tEnd
runSimulation(tStart, tEnd, nFrames, initDt)

-- print some timing information
log(string.format("Total time in vlasov solver for electrons = %g", 
		  vlasovSolverElc:totalAdvanceTime()))
log(string.format("Total time in vlasov solver for ions = %g",
		  vlasovSolverIon:totalAdvanceTime()))
log(string.format("Total time EM solver = %g",
		  maxwellSlvr:totalAdvanceTime()))
log(string.format("Total time momentum computations (elc+ion) = %g", 
		  momentumCalcElc:totalAdvanceTime()+
		     momentumCalcIon:totalAdvanceTime()))