-- Gkyl --------------------------------------------------------------
-- Collisional relaxation test - two drifting Maxwellians ------------
-- Changeset c9ca415 -------------------------------------------------
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
-- Normalization in two velocity dimensions is 1/(2*math.pi*vth^2).
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 5.0, -- End time.
   nFrame = 100, -- Number of frames to write.
   lower = {0.0}, -- Configuration space lower left.
   upper = {1.0}, -- Configuration space upper right.
   cells = {1}, -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder = 2, -- Polynomial order.
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
      lower = {-8.0, -8.0},
      upper = {8.0, 8.0},
      cells = {16, 16},
      -- Initial conditions.
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(0.5, vx, vy, 3.0, 0.0, 0.5) +
	    maxwellian2D(0.5, vx, vy, 0.0, 3.0, 0.5)
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
      diagnosticMoments = { "M0", "M1i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
}
-- Run application.
sim:run()
