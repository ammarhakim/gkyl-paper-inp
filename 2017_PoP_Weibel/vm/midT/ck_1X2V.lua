-- Gkeyll
-- Continuum kinetic 1X2V simulation
-- Current filamentation instability (Weibel)
-- Petr

----------------------------------------------------------------------
-- Simulation parameters ---------------------------------------------
pi = Lucee.Pi
-- Constatnts
chargeElc = -1.0
chargeIon = 1.0
massElc = 1.0
massIon = 1836.0
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
nIon0 = 1.0
uxIon0 = 0.0
uyIon0 = 0.0
uzIon0 = 0.0
TIon0 = 0.01
k0 = 0.4
perturb = 1e-3
-- IC automatically calculated
vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)
vthIon0 = math.sqrt(TIon0/massIon)

-- Domain and time
polyOrder = 2
cuts = {16, 2, 2}
numCellsConf = {128}
numCellsPhase = {numCellsConf[1], 48, 48} -- Nx, Nvx, Nvy
lowerBndConf = {0.0}
upperBndConf = {2*pi/k0}
lowerBndPhaseElc = {lowerBndConf[1], -1.0, -1.0} -- x, vx, vy
upperBndPhaseElc = {upperBndConf[1], 1.0, 1.0} -- x, vx, vy
lowerBndPhaseIon = {lowerBndConf[1], -6*vthIon0, -6*vthIon0} -- x, vx, vy
upperBndPhaseIon = {upperBndConf[1], 6*vthIon0, 6*vthIon0} -- x, vx, vy
periodicDirs = {0}
cfl = (1.0/3.0)/(2*polyOrder + 1)
tStart = 0.0
tEnd = 150.0
numFrames = 300

-- Other setup
Lucee.IsRestarting = false
Lucee.RestartFrame = -1
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
cellSize = (upperBndConf[1]-lowerBndConf[1])/numCellsConf[1]
log("Setup verification:")
log(" * tEnd = %g,  numFrames = %d", tEnd, numFrames)
log(" * Speed of light = %g", lightSpeed)
log(" * Time-step from light speed = %g", cfl*cellSize/lightSpeed)
log(" * Vthe1/c = %g", vthElc10/lightSpeed)
log(" * Vthe2/c = %g", vthElc20/lightSpeed)
log(" * Cell size = %g", cellSize)


----------------------------------------------------------------------
-- Computational domain, Data fields ---------------------------------
phaseDecomp = DecompRegionCalc3D.CartProd {cuts = cuts}
confDecomp = DecompRegionCalc1D.SubCartProd3D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

phaseGridElc = Grid.RectCart3D {
   lower = lowerBndPhaseElc,
   upper = upperBndPhaseElc,
   cells = numCellsPhase,
   decomposition = phaseDecomp,
   periodicDirs = periodicDirs,   
}
phaseGridIon = Grid.RectCart3D {
   lower = lowerBndPhaseIon,
   upper = upperBndPhaseIon,
   cells = numCellsPhase,
   decomposition = phaseDecomp,
   periodicDirs = periodicDirs,
}
confGrid = Grid.RectCart1D {
   lower = lowerBndConf,
   upper = upperBndConf,
   cells = numCellsConf,
   decomposition = confDecomp,
   periodicDirs = periodicDirs,
}

phaseBasisElc = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridElc,
   polyOrder = polyOrder,
}
phaseBasisIon = NodalFiniteElement3D.SerendipityElement {
   onGrid = phaseGridIon,
   polyOrder = polyOrder,
}
confBasis = NodalFiniteElement1D.LagrangeTensor {
   onGrid = confGrid,
   polyOrder = polyOrder,
   nodeLocation = "uniform",
}

createFieldElc = function()
   return DataStruct.Field3D {
      onGrid = phaseGridElc,
      numComponents = phaseBasisElc:numNodes(),
      ghost = {1, 1},
   }
end
distfElc = createFieldElc()
distfNewElc = createFieldElc()
distf1Elc = createFieldElc()
distfElcDup = createFieldElc()

createFieldIon = function()
   return DataStruct.Field3D {
      onGrid = phaseGridIon,
      numComponents = phaseBasisIon:numNodes(),
      ghost = {1, 1},
   }
end
distfIon = createFieldIon()
distfNewIon = createFieldIon()
distf1Ion = createFieldIon()

createFieldConf = function(numComponents)
   return DataStruct.Field1D {
      onGrid = confGrid,
      numComponents = numComponents*confBasis:numNodes(),
      ghost = {1, 1},
   }
end
numDensityElc = createFieldConf(1)
momentumElc = createFieldConf(2)
ptclEnergyElc = createFieldConf(3)
scalarPtclEnergyElc = createFieldConf(1)

numDensityIon = createFieldConf(1)
momentumIon = createFieldConf(2)

current = createFieldConf(2)
emSource = createFieldConf(8)
em = createFieldConf(8)
emNew = createFieldConf(8)
em1 = createFieldConf(8)
emDup = createFieldConf(8)

--emNox = createFieldConf(8)

----------------------------------------------------------------------
-- Boundary conditions -----------------------------------------------
function applyDistfBc(fld, tCurr, tAdv)
   fld:sync()
end

function applyEmBc(fld, tCurr, tAdv)
   fld:sync()
end

----------------------------------------------------------------------
-- Initial conditions setup ------------------------------------------
function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*pi*vth^2)*math.exp(-v2/(2*vth^2))
end

initDistfElc = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridElc,
   basis = phaseBasisElc,
   shareCommonNodes = false, 
   evaluate = function(x, vx, vy, t)
      return maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +
	 maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)

   end
}
initDistfIon = Updater.ProjectOnNodalBasis3D {
   onGrid = phaseGridIon,
   basis = phaseBasisIon,
   shareCommonNodes = false,
   evaluate = function(x, vx, vy, t)
      return maxwellian2D(nIon0, vx, vy, uxIon0, uyIon0, vthIon0)
   end
}
initField = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false,
   evaluate = function (x, y, z, t)
      local Ex = 0.0
      local Ey = 0.0
      local Ez = 0.0
      local Bx = 0.0
      local By = 0.0
      local Bz = perturb*math.sin(k0*x)
      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
   end
}

function initNox(x, y, z)
   return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end


----------------------------------------------------------------------
-- Equation solvers --------------------------------------------------
vlasovSolverElc = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   cfl = cfl,
   charge = chargeElc,
   mass = massElc,
   polyOrder = polyOrder,
}
vlasovSolverIon = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   cfl = cfl,
   charge = chargeIon,
   mass = massIon,
   polyOrder = polyOrder,
}

maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
   numericalFlux = "upwind",  -- "upwind" or "central"
}
maxwellSlvr = Updater.NodalDgHyper1D {
   onGrid = confGrid,
   basis = confBasis,
   equation = maxwellEqn,
   cfl = cfl,
}

scalingPositivityLimiter = Updater.ScalingPositivityLimiter1X2V {
   onGrid = phaseGridElc,
   basis  = phaseBasisElc,
}


----------------------------------------------------------------------
-- Moment updaters ---------------------------------------------------
createMomentUpdaterElc = function(moment, scalarPtclEnergy)
   if not scalarPtclEnergy then
      scalarPtclEnergy = false
   end
   return Updater.DistFuncMomentCalc1X2V {
      onGrid = phaseGridElc,
      phaseBasis = phaseBasisElc,
      confBasis = confBasis,   
      moment = moment,
      scalarPtclEnergy = scalarPtclEnergy,
   }
end
numDensityCalcElc = createMomentUpdaterElc(0)
momentumCalcElc = createMomentUpdaterElc(1)
ptclEnergyCalcElc = createMomentUpdaterElc(2)
scalarPtclEnergyCalcElc = createMomentUpdaterElc(2, true)

createMomentUpdaterIon = function(moment, scalarPtclEnergy)
   if not scalarPtclEnergy then
      scalarPtclEnergy = false
   end
   return Updater.DistFuncMomentCalc1X2V {
      onGrid = phaseGridIon,
      phaseBasis = phaseBasisIon,
      confBasis = confBasis,   
      moment = moment,
      scalarPtclEnergy = scalarPtclEnergy,
   }
end
numDensityCalcIon = createMomentUpdaterIon(0)
momentumCalcIon = createMomentUpdaterIon(1)


----------------------------------------------------------------------
-- Source updaters ---------------------------------------------------
copyToEmSource = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {0, 1},
   targetComponents = {0, 1},
}
-- This strange looking updater copies the currents into the EM source
-- field. Perhaps this is not the best way to do things, and one can
-- imagine a source updater which adds current sources to the dE/dt
-- Maxwell equation

copyB = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {3, 4, 5, 6, 7},
   targetComponents = {3, 4, 5, 6, 7},
}

----------------------------------------------------------------------
-- Solver utilities --------------------------------------------------
function runUpdater(updater, tCurr, dt, inpFlds, outFlds)
   updater:setCurrTime(tCurr)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(tCurr + dt)
end

function calcMomentsElc(tCurr, dt, distfElcIn)
   runUpdater(numDensityCalcElc, tCurr, dt, {distfElcIn}, {numDensityElc})
   runUpdater(momentumCalcElc, tCurr, dt, {distfElcIn}, {momentumElc})
   runUpdater(ptclEnergyCalcElc, tCurr, dt, {distfElcIn}, {ptclEnergyElc})
end


----------------------------------------------------------------------
-- Time-stepping -----------------------------------------------------
function rkStage(tCurr, dt, elcIn, ionIn, emIn, elcOut, emOut)
   --runUpdater(copyB, tCurr, dt, {emIn}, {emNox})

   local stElc, dtElc = runUpdater(vlasovSolverElc, tCurr, dt, 
				   {elcIn, emIn}, {elcOut})
   local stEm, dtEm = runUpdater(maxwellSlvr, tCurr, dt, {emIn}, {emOut})
   
   runUpdater(scalingPositivityLimiter, tCurr, dt, {}, {elcOut})

   runUpdater(momentumCalcElc, tCurr, dt, {elcIn}, {momentumElc})
   current:combine(chargeElc, momentumElc, chargeIon, momentumIon)
   runUpdater(copyToEmSource, tCurr, dt, {current}, {emSource})
   emOut:accumulate(-dt/epsilon0, emSource)
   
   if (stElc == false) or (stIon == false) or (stEm == false)  then
      return false, math.min(dtElc, dtEm)
   end
   return true, math.min(dtElc, dtEm)
end

function rk3(tCurr, dt)
   local status, dtSuggested
   -- RK stage 1
   status, dtSuggested = rkStage(tCurr, dt, distfElc, distfIon, em,
				     distf1Elc, em1)
   if (status == false)  then
      return false, dtSuggested
   end
   applyDistfBc(distf1Elc, tCurr, tCurr + dt)
   applyEmBc(em1, tCurr, tCurr + dt)

   -- RK stage 2
   status, dtSuggested = rkStage(tCurr, dt, distf1Elc, distf1Ion, em1,
				 distfNewElc, emNew)
   if (status == false)  then
      return false, dtSuggested
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   em1:combine(3.0/4.0, em, 1.0/4.0, emNew)
   applyDistfBc(distf1Elc, tCurr, tCurr + dt)
   applyEmBc(em1, tCurr, tCurr + dt)

   -- RK stage 3
   status, dtSuggested = rkStage(tCurr, dt, distf1Elc, distf1Ion, em1,
				 distfNewElc, emNew)
   if (status == false)  then
      return false, dtSuggested
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   em1:combine(1.0/3.0, em, 2.0/3.0, emNew)
   applyDistfBc(distf1Elc, tCurr, tCurr + dt)
   applyEmBc(em1, tCurr, tCurr + dt)   

   distfElc:copy(distf1Elc)
   em:copy(em1)
   return true, dtSuggested
end


----------------------------------------------------------------------
-- Diagnostics and Output --------------------------------------------
emEnergy = DataStruct.DynVector {numComponents = 1}
emEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*(ex*ex + ey*ey + ez*ez) +
	 0.5/mu0*(bx*bx + by*by + bz*bz)
   end,
}
emEnergyCalc:setIn({em})
emEnergyCalc:setOut({emEnergy})

ExEnergy = DataStruct.DynVector {numComponents = 1}
ExEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*ex*ex
   end,
}
ExEnergyCalc:setIn({em})
ExEnergyCalc:setOut({ExEnergy})

EyEnergy = DataStruct.DynVector {numComponents = 1}
EyEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*ey*ey
   end,
}
EyEnergyCalc:setIn( {em} )
EyEnergyCalc:setOut({EyEnergy})

EzEnergy = DataStruct.DynVector {numComponents = 1}
EzEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5*epsilon0*ez*ez
   end,
}
EzEnergyCalc:setIn({em})
EzEnergyCalc:setOut({EzEnergy})

BxEnergy = DataStruct.DynVector {numComponents = 1}
BxEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5/mu0*bx*bx
   end,
}
BxEnergyCalc:setIn({em})
BxEnergyCalc:setOut({BxEnergy})

ByEnergy = DataStruct.DynVector {numComponents = 1}
ByEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5/mu0*by*by
   end,
}
ByEnergyCalc:setIn({em})
ByEnergyCalc:setOut({ByEnergy})

BzEnergy = DataStruct.DynVector {numComponents = 1}
BzEnergyCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
      return 0.5/mu0*bz*bz
   end,
}
BzEnergyCalc:setIn({em})
BzEnergyCalc:setOut({BzEnergy})

totalPtclEnergyElc = DataStruct.DynVector {numComponents = 1}
totalPtclEnergyElcCalc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   integrand = function (energy)
      return massElc*energy
   end,
}
totalPtclEnergyElcCalc:setIn({scalarPtclEnergyElc})
totalPtclEnergyElcCalc:setOut({totalPtclEnergyElc})

function calcDiagnostics(tCurr, tAdv)
   runUpdater(scalarPtclEnergyCalcElc, tCurr, tAdv-tCurr,
	      {distfElc}, {scalarPtclEnergyElc})

   for i,diag in ipairs({emEnergyCalc,
			 ExEnergyCalc, EyEnergyCalc, EzEnergyCalc,
			 BxEnergyCalc, ByEnergyCalc, BzEnergyCalc,
			 totalPtclEnergyElcCalc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tAdv)
   end
end

function writeFields(frame, tCurr)
   distfElc:write(string.format("distfElc_%03d.h5", frame), tCurr)
   em:write(string.format("em_%03d.h5", frame), tCurr)   

   calcMomentsElc(tCurr, 0.0, distfElc)
   numDensityElc:write(string.format("numDensityElc_%03d.h5", frame), tCurr)
   momentumElc:write(string.format("momentumElc_%03d.h5", frame), tCurr)
   ptclEnergyElc:write(string.format("ptclEnergyElc_%03d.h5", frame), tCurr)

   -- diagnostics
   emEnergy:write(string.format("emEnergy_%03d.h5", frame))
   ExEnergy:write(string.format("ExEnergy_%03d.h5", frame))
   EyEnergy:write(string.format("EyEnergy_%03d.h5", frame))
   EzEnergy:write(string.format("EzEnergy_%03d.h5", frame))
   BxEnergy:write(string.format("BxEnergy_%03d.h5", frame))
   ByEnergy:write(string.format("ByEnergy_%03d.h5", frame))
   BzEnergy:write(string.format("BzEnergy_%03d.h5", frame))
   totalPtclEnergyElc:write(string.format("totalPtclEnergyElc_%03d.h5",
					  frame))
end

function writeDead(tCurr)
   distfElc:write( string.format("distfElc_dead.h5"), tCurr)
end

function logStep(step, stepFrame, frame, tCurr, dt, tStart, tEnd, tStartOS)
   local percent = (tCurr - tStart)/(tEnd - tStart)
   local remaining = os.difftime(os.time(), tStartOS)/percent*(1-percent)
   local h = math.floor(remaining/3600)
   remaining = remaining%3600
   local m = math.floor(remaining/60)
   local s = remaining%60
   log("Step %d (%d in frame %d), t=%8g, dt=%8g, %d%% %dh%02dm%02ds",
       step, stepFrame, frame, tCurr, dt, percent*100, h, m, s)
end


----------------------------------------------------------------------
-- Main function -----------------------------------------------------
function runSimulation(tStart, tEnd, numFrames, initDt)
   local frame = 1
   local tFrame = (tEnd-tStart)/numFrames
   local tNextFrame = tFrame
   local step = 1
   local stepFrame = 1
   local tCurr = tStart
   local dt = initDt
   local status, dtSuggested
   local tStartOS = os.time()
   
   if Lucee.IsRestarting then
      fileName = "q_" .. Lucee.RestartFrame .. ".h5"
      tCurr = q:read(fileName)
      if not tCurr then
         tCurr = tStart + (tEnd - tStart) * Lucee.RestartFrame / numFrames
      end
      frame = Lucee.RestartFrame + 1
      tNextFrame = tCurr + tFrame
      log('\nRestarting from frame %d tCurr = %g\n',
	  Lucee.RestartFrame, tCurr)
   else
      runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
      log("Electron IC updater took = %g", 
	  initDistfElc:totalAdvanceTime())
      runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
      log("Ion IC updater took = %g",
	  initDistfIon:totalAdvanceTime())
      runUpdater(initField, 0.0, 0.0, {}, {em})

      -- ions are not evolved
      runUpdater(numDensityCalcIon, 0, 0, {distfIon}, {numDensityIon})
      numDensityIon:write("numDensityIon.h5", 0)
      runUpdater(momentumCalcIon, 0, 0, {distfIon}, {momentumIon})
      momentumIon:write("momentumIon.h5", 0)

      --emNox:set(initNox)
   end
   applyDistfBc(distfElc, tStart, tStart)
   applyEmBc(em, tStart, tStart)
   
   if not Lucee.IsRestarting then
      calcDiagnostics(0.0, 0.0)
      writeFields(0, 0.0)
   end
   
   while true do
      distfElcDup:copy(distfElc)
      emDup:copy(em)
      
      if (tCurr + dt > tEnd) then
	 dt = tEnd - tCurr
      end
      logStep(step, stepFrame, frame, tCurr, dt,
	      tStart, tEnd, tStartOS)
      status, dtSuggested = rk3(tCurr, dt)
      
      if (status == false) then
	 -- time-step too large
	 log(" ** dt %g too large! Will retake with dt %g", dt, dtSuggested)
	 dt = dtSuggested
	 distfElc:copy(distfElcDup)
	 em:copy(emDup)
      else
	 -- check if a nan occured
	 if (distfElc:hasNan()) then
	    log(" ** NaN occured at %g! Stopping simulation", tCurr + dt)
	    writeDead(tCurr + dt)
	    break
	 end
	 
	 calcDiagnostics(tCurr, tCurr + dt)
	 tCurr = tCurr + dt
	 dt = dtSuggested
	 step = step + 1
	 stepFrame = stepFrame + 1
	 
	 if (tCurr > tNextFrame or tCurr + dt >= tEnd) then
	    log (" Writing data at time %g (frame %d) ...\n",
		 tCurr + dt, frame)
	    writeFields(frame, tCurr + dt)
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


----------------------------------------------------------------------
-- Timming -----------------------------------------------------------
log("Total time in vlasov solver for electrons = %g", 
		  vlasovSolverElc:totalAdvanceTime())
log("Total time in vlasov solver for ions = %g",
		  vlasovSolverIon:totalAdvanceTime())
log("Total time EM solver = %g",
		  maxwellSlvr:totalAdvanceTime())
log("Total time momentum computations (elc+ion) = %g", 
		  momentumCalcElc:totalAdvanceTime()+
		     momentumCalcIon:totalAdvanceTime())
