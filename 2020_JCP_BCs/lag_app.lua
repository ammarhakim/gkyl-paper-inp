-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

nCells = 4

app = Plasma.App {
   logToFile = false,

   tEnd = 1.0, -- end time
   nFrame = 1, -- number of output frames
   lower = {0}, -- configuration space lower left
   upper = {1}, -- configuration space upper right
   cells = {nCells}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   noFix = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-7.0, -7.0},
      upper = {7.0, 7.0},
      cells = {nCells, nCells},
      decompCuts = {1, 1},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 0.5*x + 0.5
	 end,
	 driftSpeed = {1.0, 0.0},
         temperature = function (t, zn)
	    local x = zn[1]
	    return 1.0 - 0.5*x
	 end,
         exactScaleM0 = false,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   scaleFix = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-7.0, -7.0},
      upper = {7.0, 7.0},
      cells = {nCells, nCells},
      decompCuts = {1, 1},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 0.5*x + 0.5
	 end,
	 driftSpeed = {1.0, 0.0},
         temperature = function (t, zn)
	    local x = zn[1]
	    return 1.0 - 0.5*x
	 end,
         exactScaleM0 = true,
         exactLagFixM012 = false,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },

   lagFix = Plasma.VlasovMaxwell.Species {
      charge = 1.0, mass = 1.0,
      -- velocity space grid
      lower = {-7.0, -7.0},
      upper = {7.0, 7.0},
      cells = {nCells, nCells},
      decompCuts = {1, 1},
      -- initial conditions
      init = Plasma.VlasovMaxwell.MaxwellianProjection {
         density = function (t, zn)
	    local x = zn[1]
	    return 0.5*x + 0.5
	 end,
	 driftSpeed = {1.0, 0.0},
         temperature = function (t, zn)
	    local x = zn[1]
	    return 1.0 - 0.5*x
	 end,
         exactScaleM0 = false,
         exactLagFixM012 = true,
      },
      evolve = true, -- evolve species?

      diagnosticMoments = { "M0", "M1i", "M2" }
   },
}
