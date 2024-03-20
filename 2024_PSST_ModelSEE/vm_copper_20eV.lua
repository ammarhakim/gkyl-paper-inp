-- Gkyl --------------------------------------------------------------
-- 1X1V copper emitting wall sheath simulation -----------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- SI units
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q0 = 1.6021766e-19
local q_e, q_i = -q0, q0
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e17, 1.0e17
local T_e, T_i = 20.0*q0, 20.0*q0

-- plasma parameters
local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local uB = math.sqrt(T_e/m_i)
local omega_pe = math.sqrt((n_e*q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0*T_e)/(n_e*q_e^2))

-- artificially decrease the speed of light
local c = 6.0*vth_e
local mu_0 = 1.0/(epsilon_0*c^2)

-- collision frequencies
local mfp = 500*lambda_D
local nu_ee = vth_e/mfp
local nu_ei = nu_ee
local nu_ii = vth_i/mfp
local nu_ie = (m_e/m_i)*nu_ee

-- collision frequency profiles
local function nu_eeProfile(t, xn)
   local x = xn[1]
   return nu_ee/(1 + math.exp(math.abs(x)/(6*lambda_D) - 8/1.5))
end
local function nu_eiProfile(t, xn)
   local x = xn[1]
   return nu_ei/(1 + math.exp(math.abs(x)/(6*lambda_D) - 8/1.5))
end
local function nu_iiProfile(t, xn)
   local x = xn[1]
   return nu_ii/(1 + math.exp(math.abs(x)/(6*lambda_D) - 8/1.5))
end
local function nu_ieProfile(t, xn)
   local x = xn[1]
   return nu_ie/(1 + math.exp(math.abs(x)/(6*lambda_D) - 8/1.5))
end

-- length of source region
local sourceLength = 40.0*lambda_D

-- initialization function
local function maxwellian(n, ux, vth, vx)
   return n/math.sqrt(2*math.pi*vth*vth)*math.exp(-(vx - ux)^2/(2*vth*vth))
end

sim = Plasma.App {
   logToFile = true,

   tEnd = 10000/omega_pe, -- end time
   nFrame = 1000, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {128.0*lambda_D}, -- configuration space upper right
   cells = {1024}, -- configuration space cells
   writeGhost = true,

   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3s4", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {128}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-4.0*vth_e},
      upper = {4.0*vth_e},
      cells = {128},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
         local x, vx = xn[1], xn[2]
         return maxwellian(n_e, 0.0, vth_e, vx)
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.ReflectBC{},
              Plasma.VlasovEmissionBC{inSpecies = {"elc"}, spectrum="gaussian", yield="furman-pivi", spectrumFit={{E0=1.971865, tau=0.883169}}, yieldFit={{mass=m_e, charge=q_e, Ehat_ts=276.8, gammahat_ts=1.885, t1=0.66, t2=0.8, t3=0.7, t4=1.0, s=1.54}}, elastic="furman-pivi", elasticFit={P1_inf=0.02, P1_hat=0.496, E_hat=1.0e-6, W=60.86, p=1.0}} },
      diagnostics = { "M0", "M1i", "M2ij", "M3i", "intM0", "Udrift", "VtSq" },
      -- plasma source
      src = Plasma.SteadySource {
         sourceSpecies = { "ion" },
         sourceLength = sourceLength,
         profile = function (t, xn)
            local x, vx = xn[1], xn[2]
            local m = maxwellian(1.0, 0.0, vth_e, vx)
            if math.abs(x) < sourceLength then
               return 2*(sourceLength - math.abs(x))/sourceLength*m
            else
               return 0.0
            end
         end,
      },

      -- species collisions
      coll = Plasma.LBOCollisions {
         collideWith = { "elc", "ion" },
         frequencies = { nu_ee, nu_ei },
      },
   },

   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-3.0*uB},
      upper = {3.0*uB},
      cells = {128},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
         local x, vx = xn[1], xn[2]
	 return maxwellian(n_i, 0.0, vth_i, vx)
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.ReflectBC{},
              Plasma.AbsorbBC{} },
      diagnostics = { "M0", "M1i", "M2ij", "M3i", "intM0", "Udrift", "VtSq" },
      -- plasma source
      src = Plasma.SteadySource {
         sourceSpecies = { "ion" },
         sourceLength = sourceLength,
         profile = function (t, xn)
            local x, vx = xn[1], xn[2]
            local m = maxwellian(1.0, 0.0, vth_i, vx)
            if math.abs(x) < sourceLength then
               return 2*(sourceLength - math.abs(x))/sourceLength*m
            else
               return 0.0
            end
         end,
      },
      -- species collisions
      coll = Plasma.LBOCollisions {
         collideWith = { "ion", "elc" },
         frequencies = { nu_ii, nu_ie },
      },
      
   },
   -- Field solver
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         local x =  xn[1]
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.Field.bcSymmetry,
              Plasma.Field.bcReflect },
   },
}
-- run application
sim:run()
