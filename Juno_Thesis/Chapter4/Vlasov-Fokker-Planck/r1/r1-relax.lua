-- Gkyl --------------------------------------------------------------
-- Collisional relaxation test - step function, polyOrder = 1 --------
-- Changeset c9ca415 -------------------------------------------------
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

sim = Plasma.App {
   logToFile = false,

   tEnd = 5.0, -- End time.
   nFrame = 100, -- Number of frames to write.
   lower = {0.0}, -- Configuration space lower left.
   upper = {1.0}, -- Configuration space upper right.
   cells = {1}, -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder = 1, -- Polynomial order.
   timeStepper = "rk3", -- One of "rk2", "rk3" or "rk3s4."
   cflFrac = 0.1, -- Additional CFL safety factor.

   -- Parallel decomposition for configuration space.
   decompCuts = {1}, -- Cuts in each configuration direction.
   useShared = false, -- Set to true to use shared memory.

   -- Specify which directions have periodic boundary conditions for configuration space.
   periodicDirs = {1},

   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0},
      upper = {6.0},
      cells = {16},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local vStep = 1.5
	 if math.abs(v) < vStep then
	    return 1.0/(2*vStep)
	 else
	    return 1.0e-6
	 end
      end,
      -- Only evolve collisions, do not evolve collisionless component.
      evolveCollisionless = false,
      evolveCollisions = true,
      -- Fokker-Planck collisions.
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut"},
	 frequencies = {1.0},
      },
      -- Diagnostics.
      diagnosticIntegratedMoments = { "intM2Flow", "intM2Thermal" },
   },
}
-- Run application.
sim:run()