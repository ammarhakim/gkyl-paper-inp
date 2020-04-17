-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").VlasovMaxwell

-- SI units
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q_e, q_i = -1.6021766e-19, 1.6021766e-19
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e17, 1.0e17
local u_e, u_i = 0.0, 0.0
local T_e, T_i = 10*1.6021766e-19, 1*1.6021766e-19

local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local uB = math.sqrt(T_e/m_i)
-- artificially decrease the speed of light
mu_0 = 1.0/(epsilon_0 * (10*vth_e)^2)
--epsilon_0 = 1.0/(mu_0 * (10*vth_e)^2)

local omega_pe = math.sqrt((n_e * q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0 * T_e)/(n_e * q_e^2))

-- initialization function
local function maxwellian(n, u, vth, v)
   return n / math.sqrt(2*math.pi*vth*vth) * 
      math.exp(-(v-u)^2/(2*vth*vth))
end

sim = Plasma.App {
   logToFile = false,

   tEnd = 120/omega_pe, -- end time
   nFrame = 12, -- number of output frames
   lower = {-128.0*lambda_D}, -- configuration space lower left
   upper = {128.0*lambda_D}, -- configuration space upper right
   cells = {256}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   cflFrac = 1.0, -- CFL "fraction". Usually 1.0
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {8}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-6.0*vth_e},
      upper = {6.0*vth_e},
      cells = {256},
      decompCuts = {1},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
	 density = n_e,
	 driftSpeed = {u_e},
	 temperature = T_e,
      },
      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcAbsorb,
              Plasma.Species.bcAbsorb },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },

   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-6.0*uB},
      upper = {6.0*uB},
      cells = {256},
      decompCuts = {1},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
	 density = n_i,
	 driftSpeed = {u_i},
	 temperature = T_i,
      },
      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcAbsorb,
              Plasma.Species.bcAbsorb },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },
   
   -- field solver
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.Field.bcReflect,
              Plasma.Field.bcReflect },
   },
}
-- run application
sim:run()
