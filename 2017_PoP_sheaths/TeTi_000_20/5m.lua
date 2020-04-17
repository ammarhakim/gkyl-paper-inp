-- 5M sheath with Poisson field solves

log = Lucee.logInfo

epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0/100. -- permiability of free space

gasGamma = 5./3
Te_Ti = 0.2 -- ratio of electron to ion temperaute
n0 = 1.0 -- initial number density
elcTemp = 1.0 -- electron temperature
elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge
B0 = 0.0 -- magnetic field

ionTemp = elcTemp/Te_Ti -- ion temperature
ionMass = 1836.2 -- ion mass
ionCharge = 1.0 -- ion charge

nNeut = 1.0 -- number density of neutral species

-- thermal speeds
ub = math.sqrt((elcTemp+ionTemp)/ionMass)
vtElc = math.sqrt(elcTemp/elcMass)
vtIon = math.sqrt(ionTemp/ionMass)
-- plasma frequency and Debye length
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
lambdaD = vtElc/wpe

-- electron and ion drift speeds
elcDrift = 0
ionDrift = 0

-- error correction speeds
elcErrorSpeedFactor = 1.0
mgnErrorSpeedFactor = 0.0

-- domain size and simulation time
XL, XU = -128*lambdaD, 128*lambdaD
tStart = 0.0 -- start time 
tEnd = 400/wpe
nFrames = 40

-- Resolution, time-stepping etc.
NX = 512
cfl =  0.9

-- cells spacing
dx = (XU-XL)/NX
dtLightSpeed = cfl*dx*math.sqrt(epsilon0*mu0)

-- ionization parameters
sigmaVr = 2*(0.5*Lucee.Pi-1)*ub/(XU-XL)
ionizationConst = sigmaVr*nNeut

-- Load number densities calculated from the Roberson model
loadedNumDensityElc = { }
loadedNumDensityIon = { }
loadedMomentumElc = { }
loadedMomentumIon = { }

local initFile = io.open("../initNumDensityElc.txt")
idx = 0
for line in initFile:lines() do
   loadedNumDensityElc[idx] = line
   loadedNumDensityElc[-idx] = line
   idx = idx+1
end
initFile:close()
local initFile = io.open("../initNumDensityIon.txt")
idx = 0
for line in initFile:lines() do
   loadedNumDensityIon[idx] = line
   loadedNumDensityIon[-idx] = line
   idx = idx+1
end
initFile:close()

local initFile = io.open("../initMomentumElc.txt")
idx = 0
for line in initFile:lines() do
   loadedMomentumElc[idx] = line
   loadedMomentumElc[-idx] = line
   loadedMomentumElc[-idx] = -1*loadedMomentumElc[-idx]
   idx = idx+1
end
initFile:close()
local initFile = io.open("../initMomentumIon.txt")
idx = 0
for line in initFile:lines() do
   loadedMomentumIon[idx] = line
   loadedMomentumIon[-idx] = line
   loadedMomentumIon[-idx] = -1*loadedMomentumIon[-idx]
   idx = idx+1
end
initFile:close()

-- print some diagnostics
log(string.format("tEnd=%g,  nFrames = %d", tEnd, nFrames))
log(string.format("Bohm speed = %g", ub))
log(string.format("Electron thermal speed = %g", vtElc))
log(string.format("Plasma frequency = %g", wpe))
log(string.format("Debye length = %g", lambdaD))
log(string.format("Cell size = %g", (math.abs(XU)+math.abs(XL))/NX))
log(string.format("Ion thermal speed = %g", vtIon))
log(string.format("Electron/Ion drift speed = %g", elcDrift))

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc1D.CartGeneral {}
-- computational domain
grid = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
   decomposition = decomp,
}

-- solution
q = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- static B field
staticEB = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 6,
   ghost = {2, 2},
}
   
-- aliases to various sub-systems
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-- potential: this is a linear DG field
phi = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 2,
   ghost = {2, 2},
}
-- densities and charge density
rhoeDG = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 2,
   ghost = {2, 2},
}
rhoiDG = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 2,
   ghost = {2, 2},
}
chargeDensity = DataStruct.Field1D {
   onGrid = grid,
   numComponents = 2,
   ghost = {2, 2},
}

-- FEM nodal basis for Poisson solve
basis = NodalFiniteElement1D.Lobatto {
   onGrid = grid,
   polyOrder = 1,
}

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

-----------------------
-- INITIAL CONDITION --
-----------------------
-- there is no round function in Lua... :( (this should work well for
-- both positive and negative numbers
function round(x)
   if x > 0 then
      return math.floor(x+0.5)
   else
      return math.ceil(x-0.5)
   end
end

-- initial conditions
function init(x,y,z)
   local idx = round(2*x)
   local rhoElc = elcMass*loadedNumDensityElc[idx]
   local rhoIon = ionMass*loadedNumDensityIon[idx]
   local momentumXElc = elcMass*loadedMomentumElc[idx]*math.sqrt(elcTemp)
   local momentumXIon = ionMass*loadedMomentumIon[idx]*math.sqrt(elcTemp)
   local erElc = loadedNumDensityElc[idx]*elcTemp/(gasGamma-1) + 
      0.5*momentumXElc*momentumXElc/rhoElc
   local erIon = loadedNumDensityIon[idx]*ionTemp/(gasGamma-1) + 
      0.5*momentumXIon*momentumXIon/rhoIon
   return rhoElc, momentumXElc, 0.0, 0.0, erElc, rhoIon, momentumXIon, 0.0, 0.0, erIon, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

-- static magnetic field
function initStaticEB(x,y,z)
   return 0, 0, 0, 0, B0, 0.0
end

-- ionization source term
function resetIonization(x, y, z)
   return 0.0, 0.0, 0.0, 0.0, 0.0
end

------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

bcElc = BoundaryCondition.Const {
   components = {0, 1, 2, 3, 4},
   values = {n0*elcMass*1e-3, 0, 0, 0, n0*elcTemp*1e-3/(gasGamma-1)}
}
bcIon = BoundaryCondition.Const {
   components = {5, 6, 7, 8, 9},
   values = {n0*ionMass*1e-3, 0, 0, 0, n0*ionTemp*1e-3/(gasGamma-1)}
}

bcRight = Updater.Bc1D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {
      bcElc, bcIon
   },
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "upper",
}
bcLeft = Updater.Bc1D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {
      bcElc, bcIon
   },
   -- direction to apply
   dir = 0,
   -- edge to apply on
   edge = "lower",
}
-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   for i,bc in ipairs({bcRight}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   for i,bc in ipairs({bcLeft}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end   
   --fld:applyCopyBc(0, "lower")

   -- sync ghost cells
   fld:sync()
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
elcEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
ionEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
elcEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
ionEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = 1.0,
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = 0.0,
}

-- ds solvers for regular Euler equations along X
elcFluidSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "minmod",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "minmod",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax Euler equations along X
elcLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxLaxSlvrDir0 = Updater.WavePropagation1D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- updater for source terms
sourceSlvr = Updater.ImplicitFiveMomentSrc1D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = 1.0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = true,
}

-- Lorentz force on electrons
elcLorentzForce = PointSource.LorentzForce {
   -- takes electron density, momentum and EM fields
   inpComponents = {0, 1, 2, 3, 10, 11, 12, 13, 14, 15},
   -- sets electron momentum and energy source
   outComponents = {1, 2, 3, 4},
   -- species charge and mass
   charge = elcCharge,
   mass = elcMass,
}
-- Lorentz force on ions
ionLorentzForce = PointSource.LorentzForce {
   -- takes ion density, momentum and EM fields
   inpComponents = {5, 6, 7, 8, 10, 11, 12, 13, 14, 15},
   -- sets ion momentum and energy source
   outComponents = {6, 7, 8, 9},
   -- species charge and mass
   charge = ionCharge,
   mass = ionMass,
}
ionization = PointSource.Ionization {
   -- takes plasma parameters and B field
   inpComponents = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15},
   -- outputs densities and energies
   outComponents = {0, 4, 5, 9},
   -- other inputs
   gasGamma = gasGamma,
   massElc = elcMass,
   massIon = ionMass,
   ionizationConst = ionizationConst,
   permeability = mu0,
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, t)
   -- two-fluid sources
   sourceSlvr:setIn( {staticEB} )
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end
-- updaters to solve ODEs for source-term splitting scheme
--sourceSlvr = Updater.GridOdePointIntegrator1D {
--   onGrid = grid,
   -- terms to include in integration step
--   terms = {elcLorentzForce, ionLorentzForce},
--}
ionizationSlvr = Updater.GridOdePointIntegrator1D {
   onGrid = grid,
   -- terms to include in integration step
   terms = {ionization},
}


-- updater to copy FV charge to DG charge field
rhoToDg = Updater.FiniteVolumeToLinearDG1D {
   onGrid = grid,
   component = 0, -- copy density
   extrapolateDomainBoundaryNodes = true, -- use extrapolation and not ghost cells
}
-- updater to compute gradient of first-order linear DG and put it on FV field
calcEx = Updater.GradLinearDGToFiniteVolume1D {
   onGrid = grid,
   component = 0, -- ex is first entry into EM field
}

-- updater to compute phi from charge density
phiFromChargeDensityCalc = Updater.FemPoisson1D {
   onGrid = grid,
   basis = basis,
   sourceNodesShared = false, -- charge density is discontinous
   solutionNodesShared = false, -- solution is  discontinous
   -- left boundary is wall, so fix potential to ground
   bcLeft = { T = "D", V = 0.0 },
   -- right boundary is wall, so fix potential to ground
   bcRight = { T = "D", V = 0.0 },
}

-- function to update ionization source terms
function updateIonization(qIn, tCurr, t)
   runUpdater(ionizationSlvr, tCurr, t-tCurr, {}, {qIn})
end

-- function to update source terms
--function updateSource(qIn, tCurr, t)
--   runUpdater(sourceSlvr, tCurr, t-tCurr, {}, {qIn})
--end

function updateField(tCurr, t, elcFluidIn, ionFluidIn, emOut)
   -- move FV data to DG fields
   runUpdater(rhoToDg, tCurr, t-tCurr, {elcFluidIn}, {rhoeDG})
   runUpdater(rhoToDg, tCurr, t-tCurr, {ionFluidIn}, {rhoiDG})

   -- compute charge density (its really rhoc/epsilon0)
   chargeDensity:combine(elcCharge/(epsilon0*elcMass), rhoeDG,
			 ionCharge/(epsilon0*ionMass), rhoiDG)
   -- solve Poisson equation with given charge density
   runUpdater(phiFromChargeDensityCalc, tCurr, t-tCurr,
	      {chargeDensity}, {phi})

   -- compute electric field by taking gradient: note that phi is
   -- assumed to have the opposite sign than is usual in EM
   runUpdater(calcEx, tCurr, t-tCurr, {phi}, {emOut})
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- Update electric field
   --updateField(tCurr, t, elcFluidNew, ionFluidNew, emFieldNew)

   return myStatus, myDtSuggested, useLaxSolver
   --return myStatus, math.min(myDtSuggested, dtLightSpeed), useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms and ionization
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   --updateSource(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)
   updateIonization(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms and ionization
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   --updateSource(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)
   updateIonization(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested, useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({elcLaxSlvrDir0, ionLaxSlvrDir0, maxLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   -- Update electric field
   --updateField(tCurr, t, elcFluidNew, ionFluidNew, emFieldNew)

   return myStatus, myDtSuggested
   --return myStatus, math.min(myDtSuggested, dtLightSpeed)
end

-- function to take one time-step with Lax Euler solver
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms and ionization
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   --updateSource(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)
   updateIonization(q, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms and ionization
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   --updateSource(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)
   updateIonization(qNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------

ionDensInCell = DataStruct.DynVector { numComponents = 1 }
elcDensInCell = DataStruct.DynVector { numComponents = 1 }

totalElcDensity = DataStruct.DynVector { numComponents = 1 }
totalIonDensity = DataStruct.DynVector { numComponents = 1 }

-- updaters
recDensInCell = Updater.RecordFieldInCell1D {
   onGrid = grid,
   cellIndex = {NX-1},
}

totalDensityCalc = Updater.IntegrateField1D {
   onGrid = grid,
   -- index of cell to record
   integrand = function(rho, rhou, rhov, rhow, er)
      return rho
   end
}

function calcDiagnostics(curr, dt)
   --runUpdater(recDensInCell, curr, dt, {rhoIonNew}, {ionDensInCell})
   --runUpdater(recDensInCell, curr, dt, {rhoElcNew}, {elcDensInCell})
   runUpdater(totalDensityCalc, curr, dt, {elcFluid}, {totalElcDensity})
   runUpdater(totalDensityCalc, curr, dt, {ionFluid}, {totalIonDensity})
   --totalDensityCalc:setIn( {elcFluid} )
   --totalDensityCalc:setOut( {totalElcDensity} )
   --totalDensityCalc:setIn( {ionFluid} )
   --totalDensityCalc:setOut( {totalIonDensity} )
end

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

-- write data to H5 files
function writeFields(frameNum, tCurr)
   q:write( string.format("q_%d.h5", frameNum), tCurr )
   --emField:write( string.format("emField_%d.h5", frameNum), tCurr )
   --phi:write( string.format("phi_%d.h5", frameNum), tCurr )
   rhoeDG:write( string.format("rhoeDG_%d.h5", frameNum), tCurr )
   rhoiDG:write( string.format("rhoiDG_%d.h5", frameNum), tCurr )
   chargeDensity:write( string.format("chargeDensity_%d.h5", frameNum), tCurr )
   ionDensInCell:write(string.format("ionDensInCell_%d.h5", frameNum), tCurr)
   elcDensInCell:write(string.format("elcDensInCell_%d.h5", frameNum), tCurr)   
   totalElcDensity:write(string.format("totalElcDensity_%d.h5", frameNum),
			 tCurr)
   totalIonDensity:write(string.format("totalIonDensity_%d.h5", frameNum),
			 tCurr)
end

----------------------------
-- TIME-STEPPING FUNCTION --
----------------------------
function runSimulation(tStart, tEnd, nFrames, initDt)

   local frame = 1
   local tFrame = (tEnd-tStart)/nFrames
   local nextIOt = tFrame
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested
   local useLaxSolver = false

   -- the grand loop 
   while true do
      -- copy q and qNew in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
        myDt = tEnd-tCurr
      end

      -- advance fluids and fields
      if (useLaxSolver) then
        -- call Lax solver if positivity violated
        log (string.format(" Taking step %5d at time %6g with dt %g (using Lax solvers)", step, tCurr, myDt))
        status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
        useLaxSolver = false
      else
        log (string.format(" Taking step %5d at time %6g with dt %g", 
			   step, tCurr, myDt))
        status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
        -- time-step too large
        log (string.format(" ** Time step %g too large! Will retake with dt %g",
			   myDt, dtSuggested))
        myDt = dtSuggested
        qNew:copy(qNewDup)
        q:copy(qDup)
      elseif (useLaxSolver == true) then
        -- negative density/pressure occured
        log (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
        q:copy(qDup)
        qNew:copy(qNewDup)
      else
        -- check if a nan occured
        if (qNew:hasNan()) then
           log (string.format(" ** NaN occured at %g! Stopping simulation",
			      tCurr))
           break
        end

        -- compute diagnostics
        calcDiagnostics(tCurr, myDt)
        -- copy updated solution back
        q:copy(qNew)
     
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
      end 
   end -- end of time-step loop
   
   return dtSuggested
end


----------------------------
-- RUNNING THE SIMULATION --
----------------------------
-- setup initial condition
q:set(init)
q:sync()
qNew:copy(q)

staticEB:set(initStaticEB)
staticEB:sync()

-- set input/output arrays for various solvers
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidNew} )
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidNew} )
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldNew} )

elcLaxSlvrDir0:setIn( {elcFluid} )
elcLaxSlvrDir0:setOut( {elcFluidNew} )
ionLaxSlvrDir0:setIn( {ionFluid} )
ionLaxSlvrDir0:setOut( {ionFluidNew} )
maxLaxSlvrDir0:setIn( {emField} )
maxLaxSlvrDir0:setOut( {emFieldNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

tStart = 0.0
initDt = 0.1
runSimulation(tStart, tEnd, nFrames, initDt)
