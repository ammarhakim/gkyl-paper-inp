local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local prng = require "sci.prng"
local rng = prng.mrg32k3a()

-- Physical constants
eps0 = 1.0
mu0 = 1.0
c = 1.0/math.sqrt(eps0*mu0)

-- Charge and mass (e=electrons, i=ions)
qe = -1.0
me = 1.0

-- Initial densities, temperatures and drift velocities
n = 1.0
ue = 0.1*c -- Keep this non-relativistic
vthe = 0.1*ue

k_max = 17.8
perturbation = 1e-3

-- Domain size and number of cells
Lx = (2*math.pi/k_max)
delta_x = 0.02*(2*math.pi/k_max)
Nx = math.ceil(Lx/delta_x)
Lve = 4*ue
delta_ve = 0.2*vthe
Nve = math.ceil(2*Lve/delta_ve)
tEnd = 100
nFrame = math.ceil(2*tEnd)

print("vthe", vthe)
print("delta_x", delta_x)
print("delta_ve", delta_ve)
print("Lx", Lx)
print("Lve", Lve)
print("Nx", Nx)
print("Nve", Nve)
print("Cells (e):", Nve*Nx)

-- Maxwellian velocity v, given the density, flow velocity and thermal speed.
local function maxwellian1D(v, n, u, vth)
   return (n/(math.sqrt(2*math.pi)*vth))*math.exp(-0.5*((v-u)/vth)^2)
end

vlasovApp = Vlasov.App {
   --------------------------------------------------------------------------------
   -- Common
   --------------------------------------------------------------------------------
   logToFile = true,

   tEnd = tEnd,           -- End time
   nFrame = nFrame,       -- Number of output frames
   lower = {0.0},         -- Lower boundary of configuration space
   upper = {Lx},          -- Upper boundary of configuration space
   cells = {Nx},          -- Configuration space cells
   basis = "serendipity", -- One of "serendipity", "maximal-order", or "tensor"
   polyOrder = 2,         -- Polynomial order
   timeStepper = "rk3s4", -- One of "rk2", "rk3", or "rk3s4"

   -- MPI decomposition for configuration space
   decompCuts = {1},      -- Cuts in each configuration direction
   useShared = true,      -- If using shared memory

   -- Boundary conditions for configuration space
   periodicDirs = {1},    -- periodic directions (both x and y)

   --------------------------------------------------------------------------------
   -- Electrons
   --------------------------------------------------------------------------------
   elc = Vlasov.Species {
      charge = qe, mass = me,
      -- Velocity space grid
      lower = {-Lve},
      upper = {Lve},
      cells = {Nve},
      -- Initial conditions
      init = function (t, xn)
         local x, vx = xn[1], xn[2]
         local fv1 = maxwellian1D(vx, n, ue, vthe)
         local fv2 = maxwellian1D(vx, n, 0, vthe) * (1 + perturbation*math.cos(k_max*x))
         return fv1 + fv2
      end,
      evolve = true,
      diagnosticMoments = {"u","vtSq"},
   },

   --------------------------------------------------------------------------------
   -- Field solver
   --------------------------------------------------------------------------------
   field = Vlasov.Field {
      epsilon0 = eps0, mu0 = mu0,
      -- hasMagneticField = false,
      init = function (t, xn)
         local x = xn[1]
         local E = -perturbation * math.sin(k_max*x) / k_max
         return E, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,
   },

}

--------------------------------------------------------------------------------
-- Run application
--------------------------------------------------------------------------------
vlasovApp:run()
