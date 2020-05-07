-- Gkyl --------------------------------------------------------------
-- Kinetic Sod-Shock, mean-free-path/Lx = 1.0/200.0 ------------------
-- Boundary conditions are periodic, and sonic point in rarefaction --
-- Changeset c9ca415 -------------------------------------------------
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- Left/right state for shock.
nl, ul, pl = 1.0, 0.75, 1.0
nr, ur, pr = 0.125, 0.0, 0.1

Lx = 1.0 -- Domain size.
mfp = Lx/200 -- Mean-free path.

-- Thermal velocity for left and right states based on left/right pressure and density.
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- Use left state as reference.
nu = vThermal/mfp -- Collision frequency.

VL, VU = -6.0*vThermal, 6.0*vThermal -- Velocity space extents.

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*math.pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 0.1, -- End time.
   nFrame = 1, -- Number of frames to write.
   lower = {-Lx}, -- Configuration space lower left.
   upper = {Lx}, -- Configuration space upper right.
   cells = {128}, -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder = 2, -- Polynomial order.
   timeStepper = "rk3", -- One of "rk2", "rk3" or "rk3s4."

   -- Parallel decomposition for configuration space.
   decompCuts = {1}, -- Cuts in each configuration direction.
   useShared = false, -- Set to true to use shared memory.

   -- Specify which directions have periodic boundary conditions for configuration space.
   periodicDirs = {1},

   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0*vThermal},
      upper = {6.0*vThermal},
      cells = {16},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local n, u, vt = nr, ur, vThermal_r
	 if math.abs(x)<0.3 then
	    n, u, vt = nl, ul, vThermal_l
	 end
	 return maxwellian(n, u, vt, v)
      end,
      -- Fokker-Planck collisions.
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut"},
	 frequencies = {nu},
      },
      -- Diagnostics.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
}
-- Run application.
sim:run()