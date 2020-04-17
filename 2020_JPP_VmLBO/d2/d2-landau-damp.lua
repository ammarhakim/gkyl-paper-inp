-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- Maxwellian in 1x1v, you may not need this but it made the example file more compact
local function maxwellian1D(n, vx, ux, mass, temp)
   local v2 = (vx - ux)^2
   return n/math.sqrt(2*math.pi*temp/mass)*math.exp(-mass*v2/(2*temp))
end

-- normalization parameters, shouldn't need to adjust
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge

nuElc = 0.05 --sample electron collision frequency

n1 = 1.0 -- number density
Te1 = 1.0 -- electron temperature

-- derived parameters
vtElc1 = math.sqrt(Te1/elcMass)
-- plasma frequency and Debye length
wpe1 = math.sqrt(elcCharge^2*n1/(epsilon0*elcMass))
lambdaD1 = vtElc1/wpe1

-- parameters for perturbation
knumber = 0.5/lambdaD1
perturbation = 1e-4

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 20.0/wpe1, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 0.9,
   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions
   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,

   -- electrons
   elc = Plasma.VlasovSpecies {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-6.0*vtElc1},
      upper = {6.0*vtElc1},
      cells = {32},
      decompCuts = {1}, -- do not change, no parallelization in velocity space currently
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
         local alpha = perturbation
         local k = knumber

	 local fv = maxwellian1D(n1, v, 0.0, elcMass, Te1) 
	 return (1+alpha*math.cos(k*x))*fv
      end,
      evolve = true, -- evolve species?
      coll = Plasma.VmLBOCollisions {
	 collideWith = { "elc" },
	 frequencies = { nuElc },
      },      
      
      -- write out density, flow, total energy, and heat flux moments
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal", "intL2" },
   },

   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local alpha = perturbation
         local k = knumber
         return -alpha*math.sin(k*xn[1])/k, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
