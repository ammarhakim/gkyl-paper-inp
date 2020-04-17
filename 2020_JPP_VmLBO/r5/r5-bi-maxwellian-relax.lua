-- Gkyl --------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- Maxwellian in 2v
local function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*vth^2)*math.exp(-v2/(2*vth^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 5.0, -- end time
   nFrame = 100, -- number of frames to write
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {1}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2", "rk3" or "rk3s4"
   cflFrac = 0.1,

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons
   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-8.0, -8.0},
      upper = {8.0, 8.0},
      cells = {16, 16},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
	 local x, vx, vy = xn[1], xn[2], xn[3]
	 return maxwellian2D(0.5, vx, vy, 3.0, 0.0, 0.5) +
	    maxwellian2D(0.5, vx, vy, 0.0, 3.0, 0.5)
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
