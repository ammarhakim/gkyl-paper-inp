-- Gkyl --------------------------------------------------------------
-- Boron nitride 1X2V material wall simulation -----------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- SI units
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q0 = 1.6021766e-19
local q_e, q_i = -q0, q0
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e17, 1.0e17
local T_e, T_i = 5*q0, 0.5*q0

-- plasma parameters
local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local uB = math.sqrt(T_e/m_i)
local omega_pe = math.sqrt((n_e*q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0*T_e)/(n_e*q_e^2))

-- artificially decrease the speed of light
local c = 6.0*vth_e
local mu_0 = 1.0/(epsilon_0*c^2)

-- collision frequencies
local mfp = 50*lambda_D
local nu_ee = vth_e/mfp
local nu_ei = nu_ee
local nu_ii = vth_i/mfp
local nu_ie = (m_e/m_i)*nu_ee

-- collision frequency profiles
local function nu_eeProfile(t, xn)
   local x = xn[1]
   return nu_ee/(1 + math.exp(x/(12*lambda_D) - 8/1.5))
end
local function nu_eiProfile(t, xn)
   local x = xn[1]
   return nu_ei/(1 + math.exp(x/(12*lambda_D) - 8/1.5))
end
local function nu_iiProfile(t, xn)
   local x = xn[1]
   return nu_ii/(1 + math.exp(x/(12*lambda_D) - 8/1.5))
end
local function nu_ieProfile(t, xn)
   local x = xn[1]
   return nu_ie/(1 + math.exp(x/(12*lambda_D) - 8/1.5))
end

-- length of source region
local sourceLength = 100.0*lambda_D

-- initialization function
local function maxwellian(n, ux, uy, vth, vx, vy)
   return n/(2*math.pi*vth*vth)*math.exp(-((vx - ux)^2 + (vy - uy)^2)/(2*vth*vth))
end

-- load data files for Roberston sheath initialization
local fh = io.open("rob_file/robertson_ne.dat")
local n_e0, i = {}, 1
for l in fh:lines() do
   n_e0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("rob_file/robertson_ni.dat")
local n_i0, i = {}, 1
for l in fh:lines() do
   n_i0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("rob_file/robertson_u.dat")
local u_0, i = {}, 1
for l in fh:lines() do
   u_0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("rob_file/robertson_E.dat")
local E_0, i = {}, 1
for l in fh:lines() do
   E_0[i] = l
   i = i + 1
end
fh.close()

sim = Plasma.App {
   logToFile = true,

   tEnd = 6000/omega_pe, -- end time
   nFrame = 600, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {1.0}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   -- nonuniform grid definition
   coordinateMap = {
      function (z) return 128.*lambda_D*(1.-(1.-z)^1.8) end
   },

   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   --cflFrac = 1.0, -- CFL "fraction". Usually 1.0
   timeStepper = "rk3s4", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {32}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-1.0, -1.0},
      upper = {1.0, 1.0},
      cells = {32, 32},
      -- nonuniform grid definition
      coordinateMap = {
         function (z) if z > 0. then return 4.*vth_e*z^1.5 else return -4.*vth_e*math.abs(z)^1.5 end end,
         function (z) if z > 0. then return 4.*vth_e*z^1.5 else return -4.*vth_e*math.abs(z)^1.5 end end
      },
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
         local x, vx, vy = xn[1], xn[2], xn[3]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
	    idx = 1280
	 else
	    idx = math.abs(x)*10/lambda_D
	    idx = math.floor(idx + 0.5) + 1
	 end
	 local n, ux = n_e0[idx]*n_e, u_0[idx]*uB*sign
         return maxwellian(n, ux, 0.0, vth_e, vx, vy)
      end,
      evolve = true, -- evolve species?
      -- vFlux = "upwind", -- use upwinding flux
      bcx = { Plasma.ReflectBC{},
              Plasma.BronoldFehskeBC{ 
                 electronAffinity = 4.5,
                 elemCharge = q0,
                 effectiveMass = 0.26,
                 electronMass = m_e 
              }
            },
      diagnostics = { "M0", "M1i", "M2ij", "M3i", "intM0", "Udrift", "VtSq" },
      -- plasma source
      src = Plasma.SteadySource {
         -- species used to calculate the outbound particle flux. The same sourceSpecies 
         -- is used for each species of particle to avoid violation of quasineutrality
         sourceSpecies = { "ion" },
         sourceLength = sourceLength, -- length of source region
         profile = function (t, xn)
            local x, vx, vy = xn[1], xn[2], xn[3]
            local m = maxwellian(1.0, 0.0, 0.0, vth_e, vx, vy)
            -- spatial profile to scale to the outbound flux (should integrate to sourceLength)
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
         frequencies = { nu_eeProfile, nu_eiProfile },
      },
   },

   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-3.0*uB, -3.0*uB},
      upper = {3.0*uB, 3.0*uB},
      cells = {32, 32},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
         local x, vx, vy = xn[1], xn[2], xn[3]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
            idx = 1280
         else
            idx = math.abs(x)*10/lambda_D
            idx = math.floor(idx + 0.5) + 1
         end
	 local n, ux = n_i0[idx]*n_i, u_0[idx]*uB*sign
         return maxwellian(n, ux, 0.0, vth_i, vx, vy)
      end,
      evolve = true, -- evolve species?
      -- vFlux = "upwind",
      bcx = { Plasma.ReflectBC{},
              Plasma.AbsorbBC{} },
      diagnostics = { "M0", "M1i", "M2ij", "M3i", "intM0", "Udrift", "VtSq" },
      -- plasma source
      src = Plasma.SteadySource {
         sourceSpecies = { "ion" },
         sourceLength = sourceLength,
         profile = function (t, xn)
            local x, vx, vy = xn[1], xn[2], xn[3]
            local m = maxwellian(1.0, 0.0, 0.0, vth_i, vx, vy)
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
         frequencies = { nu_iiProfile, nu_ieProfile },
      }, 
   },
   
   -- field solver
   field = Plasma.Field {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         local x =  xn[1]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
            idx = 1280
         else
            idx = math.abs(x)*10/lambda_D
            idx = math.floor(idx + 0.5) + 1
         end
	 local Ex = E_0[idx]*T_e/(q0*lambda_D)*sign
         return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.Field.bcSymmetry,
              Plasma.Field.bcReflect },
   },
}
-- run application
sim:run()
