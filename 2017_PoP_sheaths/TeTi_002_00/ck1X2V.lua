-- Vlasov-Maxwell solver for 1X1V sheath

----------------------------------
-- Problem dependent parameters --
----------------------------------

log = Lucee.logInfo

polyOrder = 2 -- polynomial order
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0/100 -- permeability of free space
lightSpeed = 1.0/math.sqrt(mu0*epsilon0) -- speed of light

B0 = 0

Te_Ti = 2.0 -- ratio of electron to ion temperaute
n0 = 1.0 -- initial number density
elcTemp = 1.0 -- electron temperature
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

ionTemp = elcTemp/Te_Ti -- ion temperature
ionMass = 1836.2 -- ion mass
ionCharge = 1.0 -- ion charge

nNeut = 1.0 -- number density of neutrals

-- thermal speeds
uB = math.sqrt((elcTemp+ionTemp)/ionMass)
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
vtNeut = vtIon
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- drift speeds
elcDrift = 0
ionDrift = 0
neutDrift = 0

-- domain size and simulation time
XL, XU = -128*lambdaD, 128*lambdaD
tStart = 0.0 -- start time 
tEnd = 400/wpe
nFrames = 40

-- Resolution, time-stepping etc.
NX = 256
NV = 16
polyOrder = 2

-- Collision operator variables
mfp = 10*lambdaD
nuElc = vtElc/mfp
nuIon = vtIon/mfp

-- Ionization variables
nuIonization = 2*(0.5*Lucee.Pi-1)*uB/(XU-XL)*nNeut

cfl =  0.5/(2*polyOrder+1)

-- compute max thermal speed to set velocity space extents
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc
-- (sound speeds are used for ions as we know that the ions should be
-- approximately supersonic at outlet)
VL_ION, VU_ION = -5.0*uB, 5.0*uB

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames = %d", tEnd, nFrames))
log(string.format("Bohm speed = %g", uB))
log(string.format("Electron thermal speed = %g", vtElc))
log(string.format("Plasma frequency = %g", wpe))
log(string.format("Debye length = %g", lambdaD))
log(string.format("Cell size = %g", (math.abs(XU)+math.abs(XL))/NX))
log(string.format("Ion thermal speed = %g", vtIon))
log(string.format("Electron drift speed = %g", elcDrift))
log(string.format("Ion drift speed = %g", ionDrift))
log(string.format("Configuration domain extents = [%g,%g]", XL, XU))
log(string.format("Electron domain extents = [%g,%g]", VL_ELC, VU_ELC))
log(string.format("Ion domain extents = [%g,%g]", VL_ION, VU_ION))
log(string.format("Electron collision frequency = %g", nuElc))
log(string.format("Ion collision frequency = %g", nuIon))

-- Load number densities calculated from the Roberson model
loadedNumDensityElc = { }
loadedNumDensityIon = { }
loadedFluxElc = { }
loadedFluxIon = { }

local initFile = io.open("initNumDensityElc.dat")
idx = 0
for line in initFile:lines() do
   loadedNumDensityElc[idx] = line
   idx = idx+1
end
initFile:close()
local initFile = io.open("initNumDensityIon.dat")
idx = 0
for line in initFile:lines() do
   loadedNumDensityIon[idx] = line
   idx = idx+1
end
initFile:close()

local initFile = io.open("initFluxElc.dat")
idx = 0
for line in initFile:lines() do
   loadedFluxElc[idx] = line
   idx = idx+1
end
initFile:close()
local initFile = io.open("initFluxIon.dat")
idx = 0
for line in initFile:lines() do
   loadedFluxIon[idx] = line
   idx = idx+1
end
initFile:close()

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
phaseDecomp = DecompRegionCalc3D.CartProd { cuts = {64,1,1} }
confDecomp = DecompRegionCalc1D.SubCartProd3D {
   decomposition = phaseDecomp,
   collectDirections = {0},
}

-- phase space grid for electrons
phaseGridElc = Grid.RectCart3D {
   lower = {XL, VL_ELC, VL_ELC},
   upper = {XU, VU_ELC, VU_ELC},
   cells = {NX, NV, NV},
   decomposition = phaseDecomp,  
}
-- phase space grid for ions
phaseGridIon = Grid.RectCart3D {
   lower = {XL, VL_ION, VL_ION},
   upper = {XU, VU_ION, VU_ION},
   cells = {NX, NV, NV},
   decomposition = phaseDecomp,
}

-- configuration space grid (same for electrons and ions)
confGrid = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
   decomposition = confDecomp, 
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
   nodeLocation = "lobatto"
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
distfMaxwellElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfBgkElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfNewElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfTempElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfMaxwellIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}
distfBgkIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}
distfTempIon = DataStruct.Field3D {
   onGrid = phaseGridIon,
   numComponents = phaseBasisIon:numNodes(),
   ghost = {1, 1},
}

-- ionization term distribution functions
distfIonizationElc = DataStruct.Field3D {
   onGrid = phaseGridElc,
   numComponents = phaseBasisElc:numNodes(),
   ghost = {1, 1},
}
distfIonizationIon = DataStruct.Field3D {
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
-- Electron flux
fluxElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion flux
fluxIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
diagPtclEnergyElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
ptclEnergyElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 3*confBasis:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
diagPtclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 2*confBasis:numNodes(),
   ghost = {1, 1},
}
ptclEnergyIon = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = 3*confBasis:numNodes(),
   ghost = {1, 1},
}

collFreqElc = DataStruct.Field1D {
   onGrid = confGrid,
   numComponents = confBasis:numNodes(),
   ghost = {1, 1},
}
collFreqIon = DataStruct.Field1D {
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
   --writeGhost = {1,1},
}
-- for RK time-stepping
emTemp = DataStruct.Field1D {
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
-- there is no round function in Lua... :( (this should work well for
-- both positive and negative numbers
function round(x)
   if x > 0 then
      return math.floor(x+0.5)
   else
      return math.ceil(x-0.5)
   end
end
-- updaters to apply initial conditions for number density
initProjectNumDensityElc = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      return loadedNumDensityElc[idx]
   end
}
initProjectNumDensityIon = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      return loadedNumDensityIon[idx]
   end
}
---- updaters to apply initial conditions for momenta
initProjectFluxElc = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      if x>= 0 then
	 return loadedFluxElc[idx], 0
      else
	 return -1*loadedFluxElc[idx], 0
      end
   end
}
initProjectFluxIon = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      if x >= 0 then
	 return loadedFluxIon[idx], 0
      else
	 return -1*loadedFluxIon[idx], 0
      end
   end
}
-- updaters to apply initial conditions for particle energies
initProjectPtclEnergyElc = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      return vtElc^2*loadedNumDensityElc[idx] + 
	 loadedFluxElc[idx]*loadedFluxElc[idx]/loadedNumDensityElc[idx], vtElc^2*loadedNumDensityElc[idx]
   end
}
initProjectPtclEnergyIon = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   -- basis functions to use
   basis = confBasis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local idx = round(2*math.abs(x))
      return vtIon^2*loadedNumDensityIon[idx] + 
        loadedFluxIon[idx]*loadedFluxIon[idx]/loadedNumDensityIon[idx], vtIon^2*loadedNumDensityIon[idx]
   end
}

-- updaters to initialize Maxwellian distribution from moments
initMaxwellianElc = Updater.MaxwellDistInit1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
}
initMaxwellianIon = Updater.MaxwellDistInit1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
}

-- BGK RHS updaters
bgkCollElc = Updater.BGKCollUpdater1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   mass = elcMass,
   elemCharge = ionCharge, -- expects + sign
   permitivity = epsilon0,
}
bgkCollIon = Updater.BGKCollUpdater1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   mass = ionMass,
   elemCharge = ionCharge,
   permitivity = epsilon0,
}

-- updaters to initialize Maxwellian distribution ro the ionization
-- operator
initIonizationDistfElc = Updater.MaxwellDistInit1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   zeroDriftOutput = true,
}
initIonizationDistfIon = Updater.MaxwellDistInit1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   zeroDriftOutput = true,
   arbitraryDensity = true,
}

-- updater to initialize EM fields
initField = Updater.ProjectOnNodalBasis1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      -- no fields initially
      return 0.0, 0.0, 0.0, 0.0, 0.0, B0, 0.0, 0.0
   end
}

----------------------
-- EQUATION SOLVERS --
----------------------
-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   cfl = cfl,
   charge = elcCharge,
   mass = elcMass,
   polyOrder = polyOrder,
}
vlasovSolverIon = Updater.EigenNodalVlasov1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,   
   cfl = cfl,
   charge = ionCharge,
   mass = ionMass,
   polyOrder = polyOrder,
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
   moment = 0, -- moment to compute
}
numDensityCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 0, --moment to compute
}

-- Updater to compute electron flux
fluxCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 1,
}
fluxCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 1,
}

-- Updater to compute energy
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridElc,
   phaseBasis = phaseBasisElc,
   confBasis = confBasis,
   moment = 2,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1X2V {
   onGrid = phaseGridIon,
   phaseBasis = phaseBasisIon,
   confBasis = confBasis,
   moment = 2,
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

copyPtclEnergy = Updater.CopyNodalFields1D {
   onGrid = confGrid,
   sourceBasis = confBasis,
   targetBasis = confBasis,
   sourceComponents = {0, 2},
   targetComponents = {0, 1},
}

-------------------------
-- Boundary Conditions --
-------------------------
-- boundary applicator objects for fluids and fields

function getRepTbl(polyOrder, val)
   if polyOrder == 1 then
      return {val, val, val, val, val, val, val, val}
   elseif polyOrder == 2 then
      return {val, val, val, val, val, val, val, val, val, val,
	      val, val, val, val, val, val, val, val, val, val}
   end
end
function getCountTbl(polyOrder)
   if polyOrder == 1 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   elseif polyOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	      10, 11, 12, 13, 14, 15, 16, 17, 18, 19}
   end
end

bcDistfZero = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}

-- BC updaters
setDistfElcBCLower = Updater.Bc3D {
   onGrid = phaseGridElc,
   --boundaryConditions = {bcFunctionElc},
   boundaryConditions = {bcDistfZero},
   dir = 0,
   edge = "lower",
}
setDistfElcBCUpper = Updater.Bc3D {
   onGrid = phaseGridElc,
   boundaryConditions = {bcDistfZero},
   dir = 0,
   edge = "upper",
}
setDistfIonBCLower = Updater.Bc3D {
   onGrid = phaseGridIon,
   --boundaryConditions = {bcFunctionElc},
   boundaryConditions = {bcDistfZero},
   dir = 0,
   edge = "lower",
}
setDistfIonBCUpper = Updater.Bc3D {
   onGrid = phaseGridIon,
   boundaryConditions = {bcDistfZero},
   dir = 0,
   edge = "upper",
}

-- apply boundary conditions
function applyBcDistf(curr, dt, elcIn, ionIn)
   -- apply BCs on domain boundaries
   runUpdater(setDistfElcBCUpper, curr, dt, {}, {elcIn})
   runUpdater(setDistfIonBCUpper, curr, dt, {}, {ionIn})
   runUpdater(setDistfElcBCLower, curr, dt, {}, {elcIn})
   runUpdater(setDistfIonBCLower, curr, dt, {}, {ionIn})
   -- sync the distribution function across processors
   elcIn:sync()
   ionIn:sync()
end

bcElcFld = BoundaryCondition.NodalDgZeroTangent1D { 
   components = {0, 1, 2},
   basis = confBasis,
}
bcMgnFld = BoundaryCondition.NodalDgZeroNormal1D { 
   components = {3, 4, 5}, 
   basis = confBasis,
}
bcPot = BoundaryCondition.NodalDgCopy1D { 
   components = {6, 7}, 
   fact = {-1, 1},
   basis = confBasis,
}

setEmBCUpper = Updater.Bc1D {
   onGrid = confGrid,
   -- boundary conditions to apply
   boundaryConditions = {bcElcFld, bcMgnFld, bcPot},
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "upper",
}
setEmBCLower = Updater.Bc1D {
   onGrid = confGrid,
   -- boundary conditions to apply
   boundaryConditions = {bcElcFld, bcMgnFld, bcPot},
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "lower",
}

function applyBcEM(curr, dt, emIn)
   -- apply BCs on domain boundaries
   runUpdater(setEmBCUpper, curr, dt, {}, {emIn})
   runUpdater(setEmBCLower, curr, dt, {}, {emIn})
   -- sync the EM field across processors
   emIn:sync()
end
      
----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

totalPtclElc = DataStruct.DynVector { numComponents = 1, }
totalPtclIon = DataStruct.DynVector { numComponents = 1, }

-- updater compute total number of electrons in domain
totalPtclCalcElc = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- for DG fields common nodes not shared
   integrand = function (n) return n end,
}
-- updater compute total number of ions in domain
totalPtclCalcIon = Updater.IntegrateNodalField1D {
   onGrid = confGrid,
   basis = confBasis,
   shareCommonNodes = false, -- for DG fields common nodes not shared
   integrand = function (n) return n end,
}

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
-- function to calculate flux density
function calcFlux(calculator, curr, dt, distfIn, fluxOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {fluxOut})
end
-- functions to compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
   -- number density
   runUpdater(numDensityCalcElc, curr, dt, {distfElcIn}, {numDensityElc})
   runUpdater(numDensityCalcIon, curr, dt, {distfIonIn}, {numDensityIon})
   -- flux
   runUpdater(fluxCalcElc, curr, dt, {distfElcIn}, {fluxElc})
   runUpdater(fluxCalcIon, curr, dt, {distfIonIn}, {fluxIon})
   -- energy
   runUpdater(ptclEnergyCalcElc, curr, dt, {distfElcIn}, {ptclEnergyElc})
   runUpdater(ptclEnergyCalcIon, curr, dt, {distfIonIn}, {ptclEnergyIon})
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
function calcDiagnostics(curr, dt)
   calcMoments(curr, dt, distfElc, distfIon)
   runUpdater(totalPtclCalcElc, curr, dt, {numDensityElc}, {totalPtclElc})
   runUpdater(totalPtclCalcIon, curr, dt, {numDensityIon}, {totalPtclIon})
end

----------------------------
-- Time-stepping routines --
----------------------------

-- take single RK step
function rkStage(curr, dt, elcIn, ionIn, emIn,
		 elcOut, ionOut, emOut)
   -- update distribution functions and homogenous Maxwell equations
   local stElc, dtElc = updateVlasovEqn(vlasovSolverElc, curr, dt, 
					elcIn, emIn, elcOut)
   local stIon, dtIon = updateVlasovEqn(vlasovSolverIon, curr, dt, 
					ionIn, emIn, ionOut)
   local stEm, dtEm = updateMaxwellEqn(curr, dt, emIn, emOut)
   if (stElc == false) or (stIon == false) or (stEM == false) then
      return false, math.min(dtElc, dtIon, dtEm)
   end

   calcMoments(curr, dt, elcIn, ionIn)
   runUpdater(copyPtclEnergy, curr, dt, {ptclEnergyElc}, 
	      {diagPtclEnergyElc})
   runUpdater(copyPtclEnergy, curr, dt, {ptclEnergyIon}, 
	      {diagPtclEnergyIon})

   -- BGK collisions
   stBgkElc, temp = runUpdater(bgkCollElc, curr, dt,
			       {elcIn,
				numDensityElc, fluxElc, diagPtclEnergyElc},
			       {distfBgkElc, collFreqElc})
   stBgkIon, temp = runUpdater(bgkCollIon, curr, dt,
			       {ionIn,
				numDensityIon, fluxIon, diagPtclEnergyIon},
			       {distfBgkIon, collFreqIon})
   if (stBgkElc) then
      elcOut:accumulate(dt, distfBgkElc)
   end
   
   if (stBgkIon) then
      ionOut:accumulate(dt, distfBgkIon)
   end

   stMElc, temp = runUpdater(initIonizationDistfElc, curr, dt,
			     {numDensityElc, fluxElc, diagPtclEnergyElc},
			     {distfIonizationElc})
   stMIon, temp = runUpdater(initIonizationDistfIon, curr, dt,
			     {numDensityIon, fluxIon, diagPtclEnergyIon,
			      numDensityElc}, 
			     {distfIonizationIon})
   -- update distf with Poisson eq, BGK and ionization
   if (stMElc) then
      elcOut:accumulate(dt*nuIonization, distfIonizationElc)
   end
   if (stMIon) then
      ionOut:accumulate(dt*nuIonization, distfIonizationIon)
   end
   applyBcDistf(curr, dt, elcOut, ionOut)

   -- get flux
   calcFlux(fluxCalcElc, curr, dt, elcOut, fluxElc)
   calcFlux(fluxCalcIon, curr, dt, ionOut, fluxIon)
   -- get current
   current:combine(elcCharge, fluxElc, ionCharge, fluxIon)  
   -- copy into EM sources
   runUpdater(copyToEmSource, curr, dt, {current}, {emSource})
   -- add in current source to Maxwell equation output
   emOut:accumulate(-dt/epsilon0, emSource)

   applyBcEM(curr, dt, emOut)
   
   return true, math.min(dtElc, dtIon, dtEm)
end

function rk3(tCurr, myDt)
   local status, dtSuggested
   -- RK stage 1
   status, dtSuggested = rkStage(tCurr, myDt, distfElc, distfIon, em,
				 distfTempElc, distfTempIon, emTemp)
   if status == false then
      return false, dtSuggested
   end

   -- RK stage 2
   status, dtSuggested = rkStage(tCurr, myDt, distfTempElc, distfTempIon,
				 emTemp, distfNewElc, distfNewIon, emNew)

   if status == false then
      return false, dtSuggested
   end
   distfTempElc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distfTempIon:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   emTemp:combine(3.0/4.0, em, 1.0/4.0, emNew)

   applyBcDistf(tCurr, myDt, distfTempElc, distfTempIon)
   applyBcEM(tCurr, myDt, emTemp)
 
   -- RK stage 3
   status, dtSuggested = rkStage(tCurr, myDt, distfTempElc, distfTempIon,
				 emTemp, distfNewElc, distfNewIon, emNew)
   if status == false then
      return false, dtSuggested
   end
   distfTempElc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distfTempIon:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   emTemp:combine(1.0/3.0, em, 2.0/3.0, emNew)

   applyBcDistf(tCurr, myDt, distfTempElc, distfTempIon)
   applyBcEM(tCurr, myDt, emTemp)

   distfElc:copy(distfTempElc)
   distfIon:copy(distfTempIon)
   em:copy(emTemp)
 
   return true, dtSuggested
end

-- make a duplicate in case we need it
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()
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
      distfDupIon:copy(distfIon)
      emDup:copy(em)
      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
        myDt = tEnd-tCurr
      end

      -- advance particles and fields
      log (string.format(" Taking step %5d at time %.5f with dt %.6f",
			 step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)
      if (status == false) then
	 -- time-step too large
	 log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 myDt = dtSuggested
	 distfElc:copy(distfDupElc)
	 distfIon:copy(distfDupIon)
	 em:copy(emDup)
      else
	 -- compute diagnostics
	 calcDiagnostics(tCurr, myDt)
	 -- write out data
	 if (tCurr+myDt > nextIOt or tCurr+myDt >= tEnd) then
	    log (string.format(" Writing data at time %g (frame %d) ...\n",
			       tCurr+myDt, frame))
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
	 --break
      end 
   end -- end of time-step loop
   return dtSuggested
end

-- Write out data frame 'frameNum' with at specified time 'tCurr'
function writeFields(frameNum, tCurr)
   -- distribution functions
   distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
   distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)
   
   distfMaxwellElc:write(string.format("distfElcM_%d.h5", frameNum), tCurr)
   distfMaxwellIon:write(string.format("distfIonM_%d.h5", frameNum), tCurr)
   -- electric field
   em:write(string.format("em_%d.h5", frameNum), tCurr)
   -- moments
   numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)
   --chargeDensity:write(string.format("chargeDensity_%d.h5", frameNum), tCurr)
   fluxElc:write(string.format("fluxElc_%d.h5", frameNum), tCurr)
   fluxIon:write(string.format("fluxIon_%d.h5", frameNum), tCurr)
   ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)
   ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr) 
   --scalarPtclEnergyIon:write(string.format("scalarPtclEnergyIon_%d.h5", frameNum), tCurr) 
   -- diagnostics
   totalPtclElc:write(string.format("totalPtclElc_%d.h5", frameNum), tCurr)
   totalPtclIon:write(string.format("totalPtclIon_%d.h5", frameNum), tCurr)

   collFreqElc:write(string.format("collFreqElc_%d.h5", frameNum), tCurr)
   collFreqIon:write(string.format("collFreqIon_%d.h5", frameNum), tCurr)
end

----------------------------
-- RUNNING THE SIMULATION --
----------------------------
runUpdater(initProjectNumDensityElc, 0.0, 0.0, {}, {numDensityElc})
runUpdater(initProjectFluxElc, 0.0, 0.0, {}, {fluxElc})
runUpdater(initProjectPtclEnergyElc, 0.0, 0.0, {}, {diagPtclEnergyElc})
runUpdater(initProjectNumDensityIon, 0.0, 0.0, {}, {numDensityIon})
runUpdater(initProjectFluxIon, 0.0, 0.0, {}, {fluxIon})
runUpdater(initProjectPtclEnergyIon, 0.0, 0.0, {}, {diagPtclEnergyIon})
runUpdater(initMaxwellianElc, 0.0, 0.0, 
	   {numDensityElc, fluxElc, diagPtclEnergyElc}, {distfElc})
runUpdater(initMaxwellianIon, 0.0, 0.0, 
	   {numDensityIon, fluxIon, diagPtclEnergyIon}, {distfIon})



--runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
--runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})
applyBcDistf(0.0, 0.0, distfElc, distfIon)

-- initialize em field
runUpdater(initField, 0.0, 0.0, {}, {em})
applyBcEM(0, 0, em)

--compute initial diagnostics
calcDiagnostics(0.0, 0.0)

-- write out initial fields
writeFields(0, 0.0)

cntElc = 0
cntElcTot = 0
cntIon = 0
cntIonTot = 0

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
log(string.format("Total time flux computations (elc+ion) = %g", 
		  fluxCalcElc:totalAdvanceTime()+
		     fluxCalcIon:totalAdvanceTime()))