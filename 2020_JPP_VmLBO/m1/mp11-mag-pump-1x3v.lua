-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell

-- Test in which we will drive a magnetic field to reproduce magnetic pumping
-- heating thanks to the pitch angle scattering of particles.

-- Maxwellian in 1x3v, you may not need this but it made the example file more compact
local function maxwellian1D(n, vx, vy, vz, ux, uy, uz, mass, temp)
   local v2 = (vx - ux)^2+(vy - uy)^2+(vz - uz)^2
   return (n/((math.sqrt(2*math.pi*temp/mass))^3))*math.exp(-mass*v2/(2*temp))
end

-- Universal constants.
epsilon0   = 1.0                         -- permittivity of free space.
mu0        = 1.0                         -- permeability of free space.
lightSpeed = 1.0/math.sqrt(mu0*epsilon0) -- speed of light.

-- User inputs.
massRatio = 1836.153
tau       = 1.0          -- Ratio of ion to electron temperature.
vA        = 0.25          -- Alven speed.
beta      = 0.01         -- Total plasma beta.

ionMass   = massRatio    -- Ion mass in simulation.
elcMass   = 1.0          -- Electron mass in simulation.
ionCharge = 1.0          -- Ion charge in simulation.
elcCharge = -1.0         -- Electron charge in simulation.

B0  = vA        -- Driven magnetic field amplitude.

n   = 0.01                             -- Plasma density, same for ions an electrons.
Te  = beta*(B0^2)/(2.0*mu0*(1.0+tau)) -- Electron temperature.
Ti  = tau*Te                          -- Ion temperature.

-- Derived parameters.
vtElc   = math.sqrt(2.0*Te/elcMass)
vtIon   = math.sqrt(2.0*Ti/ionMass)

-- cyclotron frequency
omegaCe = ionCharge*B0/elcMass

-- gyro radius
rhoe = vtElc/omegaCe

-- Frequency of time-varying driven Bz.
w0 = 0.1*omegaCe

-- ramp time
t_ramp = 2.0*math.pi/w0

-- collision frequencies
nuElc = 0.1*w0*(math.pi^2)                     -- Sample electron collision frequency.
nuIon = nuElc/math.sqrt(ionMass*(tau^3))   -- Sample ion collision frequency.

-- Domain size.
Lx = 200.0*math.pi*rhoe
Nx = 256
-- number of cores
Nc = 128

-- points to drive current at
point1 = 0.25*Lx
point2 = 0.75*Lx
-- dx for spread of gaussian
dx = Lx/Nx

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd         = 100.0*math.pi/omegaCe,       -- End time.
   nFrame       = 400,                   -- Number of output frames.
   lower        = {0.0},               -- Configuration space lower left.
   upper        = {Lx},                -- Configuration space upper right.
   cells        = {Nx},                -- Configuration space cells.
   basis        = "serendipity",       -- One of "serendipity" or "maximal-order".
   polyOrder    = 2,                   -- Polynomial order.
   timeStepper  = "rk3",               -- One of "rk2" or "rk3".
   -- Boundary conditions for configuration space.
   periodicDirs = {1},                 -- Periodic directions.
   -- decomposition for configuration space
   decompCuts   = {Nc},                 -- Cuts in each configuration direction.
   useShared    = false,               -- If to use shared memory.

   -- Integrated moment flag, compute quantities 1000 times in simulation.
   calcIntQuantEvery = 0.001,

   -- electrons
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- Velocity space grid.
      lower      = {-8.0*vtElc,-8.0*vtElc,-8.0*vtElc},
      upper      = { 8.0*vtElc, 8.0*vtElc, 8.0*vtElc},
      cells      = {24,24,24},
      decompCuts = {1,1,1}, -- Do not change, no parallelization in velocity space currently.
      -- Initial conditions.
      init = function (t, xn)
         local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]
         local ux, uy, uz    = 0.0, 0.0, 0.0
	 local fv   = maxwellian1D(n, vx, vy, vz, ux, uy, uz, elcMass, Te) 
	 return fv
      end,
      nDistFuncFrame = 40,
      evolve = true, -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments           = { "M0", "M1i", "M2" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      coll = Plasma.LBOCollisions {
            collFreq = nuElc,
      },
   },

   -- protons
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- Velocity space grid.
      lower      = {-8.0*vtIon,-8.0*vtIon,-8.0*vtIon},
      upper      = { 8.0*vtIon, 8.0*vtIon, 8.0*vtIon},
      cells      = {24,24,24},
      decompCuts = {1,1,1}, -- Do not change, no parallelization in velocity space currently.
      -- Initial conditions.
      init = function (t, xn)
         local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]
         local ux, uy, uz    = 0.0, 0.0, 0.0
	 local fv   = maxwellian1D(n, vx, vy, vz, ux, uy, uz, ionMass, Ti)
	 return fv
      end,
      nDistFuncFrame = 40,
      evolve = true, -- Evolve species?
      -- Write out density, flow, total energy, and heat flux moments.
      diagnosticMoments           = { "M0", "M1i", "M2" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
      coll = Plasma.LBOCollisions {
            collFreq = nuIon,
      },
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
        local x = xn[1]
        return 0.0, 0.0, 0.0, 0.0, 0.0, B0
      end,
      evolve = true, -- evolve field?
   },

   -- current antenna
   driveSpecies = Plasma.FuncSpecies {
      charge = 1.0, mass = 1.0,
      momentumDensity = function (t, xn)
         local x = xn[1]
         local Jy = 0.0
	 local J0 = 0.5*B0*math.sin(0.5*math.pi*math.min(1, t/t_ramp))^2*math.sin(w0*t)
	 Jy = J0*(math.exp(-(x-point1)^2/(2.0*dx)^2)-math.exp(-(x-point2)^2/(2.0*dx)^2))
         return 0.0, Jy, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
