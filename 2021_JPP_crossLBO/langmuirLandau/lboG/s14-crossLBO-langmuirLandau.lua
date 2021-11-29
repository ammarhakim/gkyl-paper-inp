-- Gkyl ------------------------------------------------------------------------
--
-- Landau damping of a Langmuir wave in 1X3V, kinetic electrons immobile ions.
--
--

local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- Normalization parameters, shouldn't need to adjust.
epsilon0   = 1.0                       -- Permittivity of free space.
mu0        = 1.0                       -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- Speed of light.

elcMass    =  1.0       -- Electron mass.
elcCharge  = -1.0       -- Electron charge.
ionMass    =  1836.0    -- Ion mass.
ionCharge  =  1.0       -- Ion charge.

mRat       = ionMass/elcMass    -- Mass ratio, m_i/m_e.

nuElc      = 0.09                    -- Electron-electron collision frequency.
nuIon      = nuElc/math.sqrt(mRat)    -- Ion-ion collision frequency.
nuElcIon   = nuElc*math.sqrt(2)       -- Electron-ion collision frequency.
nuIonElc   = nuElcIon/mRat            -- Ion-electron collision frequency.

n0  = 1.0    -- Number density.
Te0 = 1.0    -- Electron temperature.
Ti0 = 1.0    -- Ion temperature.

-- Derived parameters
vtElc   = math.sqrt(Te0/elcMass)                            -- Electron thermal speed.
vtIon   = math.sqrt(Ti0/ionMass)                            -- Ion thermal speed.
wpe     = math.sqrt((elcCharge^2)*n0/(epsilon0*elcMass))    -- Plasma frequency.
lambdaD = vtElc/wpe                                         -- Debye length.

-- Parameters for perturbation.
kNumber       = 0.3/lambdaD
pertAmplitude = 1e-4

-- Maxwellian in 1x1v.
local function maxwellian3V(n, vx, vy, vz, ux, uy, uz, vt)
   local vSq = (vx - ux)^2+(vy - uy)^2+(vz - uz)^2
   return (n/(math.sqrt(2*math.pi*(vt^2))^3))*math.exp(-vSq/(2*(vt^2)))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd         = 100.0/wpe,              -- End time.
   nFrame       = 100,                   -- Number of output frames.
   lower        = {-math.pi/kNumber},    -- Configuration space lower left.
   upper        = { math.pi/kNumber},    -- Configuration space upper right.
   cells        = {16},                  -- Configuration space cells.
   basis        = "serendipity",         -- One of "serendipity" or "maximal-order".
   polyOrder    = 2,                     -- Polynomial order.
   cflFrac      = 0.5,
   timeStepper  = "rk3",                 -- One of "rk2" or "rk3".
   -- Boundary conditions for configuration space.
   periodicDirs = {1},                   -- Periodic directions.
   -- Decomposition for configuration space.
   decompCuts   = {8},                  -- Cuts in each configuration direction, or nodes (useShared=True).
   useShared    = true,                 -- If to use shared memory.
   restartFrameEvery = 0.02,

   -- Integrated moment flag, aim to compute 10x times than frames.
   calcIntQuantEvery = 1./1000.,

   -- Electrons.
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- Velocity space grid.
      lower = {-5.0*vtElc, -5.0*vtElc, -5.0*vtElc},
      upper = { 5.0*vtElc,  5.0*vtElc,  5.0*vtElc},
      cells = {36, 36, 36},
      vFlux = "upwind",  -- Use upwind fluxes in velocity space.
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]
         local alpha = pertAmplitude
         local k     = kNumber

	 local fv = maxwellian3V(n0, vx, vy, vz, 0.0, 0.0, 0.0, vtElc) 
	 return (1.0+alpha*math.cos(k*x))*fv
      end,
      evolve = true,    -- Evolve species?
      nDistFuncFrame = 1,
      coll = Plasma.LBOCollisions {
         collideWith = { 'elc', 'ion' },
         frequencies = { nuElc, nuElcIon },
         betaGreene  = -0.65,
      },
   },

   -- Ions.
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- Velocity space grid.
      lower = {-5.0*vtIon, -5.0*vtIon, -5.0*vtIon},
      upper = { 5.0*vtIon,  5.0*vtIon,  5.0*vtIon},
      cells = {32, 32, 32},
      vFlux = "upwind",  -- Use upwind fluxes in velocity space.
      -- Initial conditions.
      init = function (t, xn)
         local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]

         return maxwellian3V(n0, vx, vy, vz, 0.0, 0.0, 0.0, vtIon)
      end,
      evolve = false,    -- Evolve species?
      nDistFuncFrame = 1,
      coll = Plasma.LBOCollisions {
         collideWith = { 'elc', 'ion' },
         frequencies = { nuIonElc, nuIon },
         betaGreene  = -0.65,
      },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = pertAmplitude
         local k     = kNumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
