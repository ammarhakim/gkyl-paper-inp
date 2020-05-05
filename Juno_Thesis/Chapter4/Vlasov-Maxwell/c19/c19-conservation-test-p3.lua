-- Gkyl ------------------------------------------------------------------
-- Conservation test -----------------------------------------------------
-- p = 3, dx = 3 lambda_{D}, dv = 1/4 v_{th}, dt = 1/20 omega_{pe}^{-1} --
-- Changeset c9ca415 -----------------------------------------------------
-- -----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- Normalization parameters.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- Speed of light.

elcMass = 1.0 -- Electron mass.
elcCharge = -1.0 -- Electron charge.
ionMass = 1836.153 -- Proton mass.
ionCharge = 1.0 -- Proton charge.
Te_Ti = 1.0 -- Ratio of electron to ion temperature.

-- Plasma parameters.
n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 1.0/Te_Ti -- Proton temperature.

-- Derived parameters
vtElc = math.sqrt(Te/elcMass) -- Electron thermal velocity.
vtIon = math.sqrt(Ti/ionMass) -- Proton thermal velocity
-- Plasma frequency and Debye length.
wpe = math.sqrt(elcCharge^2*n0/(epsilon0*elcMass))
wpi = math.sqrt(ionCharge^2*n0/(epsilon0*ionMass))
lambdaD = vtElc/wpe

-- Full-width at half maximum (FWHM) of density profile on left and right.
-- Beta = 2 sigma^2 where sigma is FWHM.
beta_left = 2*lambdaD^2 -- FWHM on left = 1 lambdaD.
beta_right = 32*lambdaD^2 -- FWHM on right = 4 lambdaD.

-- Domain size and simulation time.
LX = 96.0*lambdaD
tEnd = 1000/wpe

sim = Plasma.App {
   logToFile = false,

   tEnd = tEnd, -- End time.
   nFrame = 1, -- Number of output frames.
   lower = {0.0}, -- Configuration space lower left.
   upper = {LX}, -- Configuration space upper right.
   cells = {32}, -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder = 3, -- Polynomial order.
   timeStepper = "rk3", -- One of "rk2", "rk3" or "rk3s4."

   -- Set size of maximum time-step.
   maximumDt = 0.05,
   -- Parallel decomposition for configuration space.
   decompCuts = {1}, -- Cuts in each configuration direction.
   useShared = false, -- Set to true to use shared memory.

   -- Specify which directions have periodic boundary conditions for configuration space.
   periodicDirs = {1},

   -- Frequency of computing integrated quantities; compute 50 times in simulation.
   calcIntQuantEvery = 0.02,

   -- Electrons.
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- Velocity space grid.
      lower = {-5.0*vtElc},
      upper = {7.0*vtElc},
      cells = {48},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vd = vtElc
	 local xm = 0.25*LX
	 -- Perturbation in the density.
	 local fv = (1.0+4.0*math.exp(-(x-xm)^2/beta_left))*maxwellian(n0, vd, vtElc, v) 
         if x>xm then
            fv = (1.0+4.0*math.exp(-(x-xm)^2/beta_right))*maxwellian(n0, vd, vtElc, v)
         end
	 return fv
      end,
      -- Diagnostics.
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },
   -- Protons.
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- Velocity space grid; centered on initial drift.
      lower = {-6.0*vtIon+vtElc},
      upper = {6.0*vtIon+vtElc},
      cells = {48},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vd = vtElc
	 local xm = 0.25*LX
	 -- Perturbation in the density.
	 local fv = (1.0+4.0*math.exp(-(x-xm)^2/beta_left))*maxwellian(n0, vd, vtIon, v) 
         if x>xm then
            fv = (1.0+4.0*math.exp(-(x-xm)^2/beta_right))*maxwellian(n0, vd, vtIon, v)
         end
	 return fv
      end,
      -- Diagnostics.
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      -- No initial electromagnetic fields.
      init = function (t, xn)
	 return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
   },
}
-- Run application.
sim:run()
