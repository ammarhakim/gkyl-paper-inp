-- Gkyl --------------------------------------------------------------
-- Kinetic Sod-Shock, mean-free-path/Lx = 1.0/500.0 ------------------
-- Changeset c9ca415 -------------------------------------------------
----------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- Left/right state for shock.
nl, ul, pl = 1.0, 0.0, 1.0
nr, ur, pr = 0.125, 0.0, 0.1

Lx = 1.0 -- Domain size.
mfp = Lx/500.0 -- Mean-free path.

-- Thermal velocity for left and right states based on left/right pressure and density.
vThermal_l = math.sqrt(pl/nl)
vThermal_r = math.sqrt(pr/nr)

vThermal = vThermal_l -- Use left state as reference.
nu = vThermal/mfp -- Collision frequency.

VL, VU = -6.0*vThermal, 6.0*vThermal

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
   lower = {0.0}, -- Configuration space lower left.
   upper = {Lx}, -- Configuration space upper right.
   cells = {64}, -- Configuration space cells.
   basis = "serendipity", -- One of "serendipity" or "maximal-order."
   polyOrder = 2, -- Polynomial order.
   timeStepper = "rk3", -- One of "rk2", "rk3" or "rk3s4."

   -- Parallel decomposition for configuration space.
   decompCuts = {1}, -- Cuts in each configuration direction.
   useShared = false, -- Set to true to use shared memory.

   -- Specify which directions have periodic oundary conditions for configuration space.
   periodicDirs = {},

   neut = Plasma.Species {
      charge = 1.0, mass = 1.0,
      -- Velocity space grid.
      lower = {-6.0*vThermal},
      upper = {6.0*vThermal},
      cells = {16},
      -- Initial conditions.
      init = function (t, xn)
	 local x, v = xn[1], xn[2]
	 local n, u, vt = nl, ul, vThermal_l
	 if x>0.5 then
	    n, u, vt = nr, ur, vThermal_r
	 end
	 return maxwellian(n, u, vt, v)
      end,
      -- Fokker-Planck collisions.
      lbo = Plasma.LBOCollisions {
	 collideWith = {"neut"},
	 frequencies = {nu},
      },
      -- Non-periodic boundary conditions.
      bcx = { Plasma.Species.bcCopy, Plasma.Species.bcCopy },
      -- Diagnostics.
      diagnosticMoments = { "M0", "M1i", "M2", "M3i" },
   },
}
-- Run application.
sim:run()