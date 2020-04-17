-- Wiebel instability

Pi = Lucee.Pi
log = Lucee.logInfo

-- physical parameters (two electron streams: "ions" are really electrons!)
gasGamma = 5./3.
elcCharge = -1.0
ionCharge = -1.0
ionMass = 1.0
elcMass = 1.0
lightSpeed = 1.0
epsilon0 = 1.0
mu0 = 1.0
mgnErrorSpeedFactor = 1.0

-- simulation parameters
udrift1 = 0.3
udrift2 = -0.3
k0 = 0.2 -- wave number
elcTemp = 0.01

Lx = 1.0
Ly = 2*Lucee.Pi/k0
n0 = 0.5

-- resolution and time-stepping
NX = 1
NY = 1024
cfl = 0.9
tStart = 0.0
tEnd = 150.0
nFrames = 75

------------------------------------------------
-- COMPUTATIONAL DOMAIN, DATA STRUCTURE, ETC. --
------------------------------------------------
-- decomposition object
decomp = DecompRegionCalc2D.CartGeneral {}
-- computational domain
grid = Grid.RectCart2D {
   lower = {0.0, 0.0},
   upper = {1.0, Ly},
   cells = {NX, NY},
   decomposition = decomp,
   periodicDirs = {0, 1},
}

-- solution
q = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- solution after update along X (ds algorithm)
qX = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- final updated solution
qNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
-- duplicate copy in case we need to take the step again
qDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}
qNewDup = DataStruct.Field2D {
   onGrid = grid,
   numComponents = 18,
   ghost = {2, 2},
}

-- aliases to various sub-systems
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidX = qX:alias(0, 5)
ionFluidX = qX:alias(5, 10)
emFieldX = qX:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function init(x,y,z)
   local rho = elcMass*n0
   local pr0 = n0*elcTemp
   local Bz = 1e-3*math.sin(k0*y)
   
   local u1 = udrift1
   local Er1 = pr0/(gasGamma-1) + 0.5*rho*u1^2

   local u2 = udrift2
   local Er2 = pr0/(gasGamma-1) + 0.5*rho*u2^2
   
   return rho, rho*u1, 0.0, 0.0, Er1, rho, rho*u2, 0.0, 0.0, Er2, 0.0, 0.0, 0.0, 0.0, 0.0, Bz, 0.0, 0.0
end


------------------------
-- Boundary Condition --
------------------------
-- boundary applicator objects for fluids and fields

-- function to apply boundary conditions to specified field
function applyBc(fld, tCurr, myDt)
   for i,bc in ipairs({}) do
      bc:setOut( {fld} )
      bc:advance(tCurr+myDt)
   end
   -- copy BCs in Y
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")   
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
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = 0.0,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-- ds solvers for regular Euler equations along X
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for regular Euler equations along Y
elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- ds solvers for Lax Euler equations along X
elcLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax Euler equations along Y
elcLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- updater for source terms
sourceSlvr = Updater.ImplicitFiveMomentSrc2D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

elcIonMomRelax = Updater.TwoFluidMomentumRelaxation2D {
   onGrid = grid,
   electronIonCollisionFrequency = 0.0,
   frictionFactor = 0.5
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, t)
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = False
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidX) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidX) == false)) then
      useLaxSolver = true
   end

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir1, ionFluidSlvrDir1, maxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   if ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
       useLaxSolver = true
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
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

   applyBc(qX, tCurr, t-tCurr)

   -- Y-direction updates
   for i,slvr in ipairs({elcLaxSlvrDir1, ionLaxSlvrDir1, maxLaxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end

   return myStatus, myDtSuggested
end

-- function to take one time-step with Lax Euler solver
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------
-- dynvector to store integrated flux

-- dynvector to store magnetic energy
bz2 = DataStruct.DynVector { numComponents = 1 }
bz2Calc = Updater.IntegrateField2D {
   onGrid = grid,
   -- index of cell to record
   integrand = function (ex, ey, ez, bx, by, bz, e1, e2)
		  return bz^2
	       end,
}
bz2Calc:setIn( {emField} )
bz2Calc:setOut( {bz2} )

-- compute diagnostic
function calcDiagnostics(tCurr, myDt)
   for i,diag in ipairs({bz2Calc}) do
      diag:setCurrTime(tCurr)
      diag:advance(tCurr+myDt)
   end
end

-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )
   bz2:write( string.format("bz2_%d.h5", frame) )
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
        log (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
        status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
        -- time-step too large
        log (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
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
           log (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
           break
        end

        -- compute diagnostics
        calcDiagnostics(tCurr, myDt)
        -- copy updated solution back
        q:copy(qNew)
     
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


----------------------------
-- RUNNING THE SIMULATION --
----------------------------
-- setup initial condition
q:set(init)
q:sync()
qNew:copy(q)

-- set input/output arrays for various solvers
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidX} )
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidX} )
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldX} )

elcFluidSlvrDir1:setIn( {elcFluidX} )
elcFluidSlvrDir1:setOut( {elcFluidNew} )
ionFluidSlvrDir1:setIn( {ionFluidX} )
ionFluidSlvrDir1:setOut( {ionFluidNew} )
maxSlvrDir1:setIn( {emFieldX} )
maxSlvrDir1:setOut( {emFieldNew} )

elcLaxSlvrDir0:setIn( {elcFluid} )
elcLaxSlvrDir0:setOut( {elcFluidX} )
ionLaxSlvrDir0:setIn( {ionFluid} )
ionLaxSlvrDir0:setOut( {ionFluidX} )
maxLaxSlvrDir0:setIn( {emField} )
maxLaxSlvrDir0:setOut( {emFieldX} )

elcLaxSlvrDir1:setIn( {elcFluidX} )
elcLaxSlvrDir1:setOut( {elcFluidNew} )
ionLaxSlvrDir1:setIn( {ionFluidX} )
ionLaxSlvrDir1:setOut( {ionFluidNew} )
maxLaxSlvrDir1:setIn( {emFieldX} )
maxLaxSlvrDir1:setOut( {emFieldNew} )

-- apply BCs on initial conditions
applyBc(q, 0.0, 0.0)
applyBc(qNew, 0.0, 0.0)

-- write initial conditions
calcDiagnostics(0.0, 0.0)
writeFields(0, 0.0)

initDt = 100.0
runSimulation(tStart, tEnd, nFrames, initDt)


