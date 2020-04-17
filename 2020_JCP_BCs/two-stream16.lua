-- Gkyl --------------------------------------------------------------
-- Electron Two-stream instability -----------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

kNumber = 0.5 -- wave-number
vth_e = 0.2 -- electron thermal velocity
vd_e = 1.0 -- drift velocity
perturbation = 1.0e-6 -- distribution function perturbation

-- initialization function
local function maxwellian1D(n, vd, vth, vx)
   return n / math.sqrt(2*math.pi*vth*vth) * 
      math.exp(-(vx-vd)^2/(2*vth*vth))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 100.0, -- end time
   nFrame = 2, -- number of output frames
   lower = {-math.pi/kNumber}, -- configuration space lower left
   upper = {math.pi/kNumber}, -- configuration space upper right
   cells = {64}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"

   -- decomposition for configuration space
   decompCuts = {2}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory
   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = -1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6},
      upper = {6},
      cells = {16},
      decompCuts = {1},
      -- initial conditions
      fp = Plasma.MaxwellianProjection {
	 density = function (t, xn)
	    local x = xn[1]
	    return 0.5*(1 + perturbation*math.cos(kNumber*x))
	 end,
	 drift = vd_e,
	 temperature = vth_e^2,
	 exactScaleM0 = false,
	 exactLagFixM012 = true,
      },
      fm = Plasma.MaxwellianProjection {
	 density = function (t, xn)
	    local x = xn[1]
	    return 0.5*(1 + perturbation*math.cos(kNumber*x))
	 end,
	 drift = -vd_e,
	 temperature = vth_e^2,
	 exactScaleM0 = false,
	 exactLagFixM012 = true,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M2" }
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x = xn[1]
	 local Ex = -perturbation * math.sin(kNumber*x) / kNumber
	 return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
   },
}
-- run application
sim:run()
