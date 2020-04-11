-- Gkyl --------------------------------------------------------------
-- Collisional relaxation test ---------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

sim = Plasma.App {
   logToFile = false,

   tEnd = 5.0, -- end time
   nFrame = 100, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {2*math.pi}, -- configuration space upper right
   cells = {2}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"
   cflFrac = 0.1/2,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-6.0},
      upper = {6.0},
      cells = {8},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vStep = 1.5
	 if math.abs(v) < vStep then
	    return 1.0/(2*vStep)
	 else
	    return 1.0e-6
	 end
      end,
      evolveCollisionless = false,
      evolveCollisions = true,
      -- collisions
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut"},
	 frequencies = {1.0},
      },
      -- diagnostics
      diagnosticMoments = { "M0", "M1i" },
      diagnosticIntegratedMoments = {
	 "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
}
-- run application
sim:run()
