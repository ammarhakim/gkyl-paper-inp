-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- SI units
local elemCharge = 1.6021766e-19
local epsilon_0, mu_0 = 8.854e-12, 1.257e-6
local q_e, q_i = -elemCharge, elemCharge
local m_e, m_i = 9.109383e-31, 1.6726218e-27
local n_e, n_i = 1.0e17, 1.0e17
local u_e, u_i = 0.0, 0.0
local T_e, T_i = 10*elemCharge, 1*elemCharge

local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local uB = math.sqrt(T_e/m_i)
local omega_pe = math.sqrt((n_e * q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0 * T_e)/(n_e * q_e^2))
-- artificially decrease the speed of light
mu_0 = 1.0/(epsilon_0 * (10*vth_e)^2)

-- initialization function
local function maxwellian(n, u, vth, v)
   return n / math.sqrt(2*math.pi*vth*vth) * 
      math.exp(-(v-u)^2/(2*vth*vth))
end

local xi = 1.0
local dv = math.sqrt((2*xi*elemCharge)/m_e)*2
local vCells = 16

-- load data files
local fh = io.open("../robertson_ne.dat")
local n_e0, i = {}, 1
for l in fh:lines() do
   n_e0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("../robertson_ni.dat")
local n_i0, i = {}, 1
for l in fh:lines() do
   n_i0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("../robertson_u.dat")
local u_0, i = {}, 1
for l in fh:lines() do
   u_0[i] = l
   i = i + 1
end
fh.close()
local fh = io.open("../robertson_E.dat")
local E_0, i = {}, 1
for l in fh:lines() do
   E_0[i] = l
   i = i + 1
end
fh.close()

sim = Plasma.App {
   logToFile = false,

   tEnd = 500/omega_pe, -- end time
   nFrame = 1, -- number of output frames
   lower = {0.0}, -- configuration space lower left
   upper = {128.0*lambda_D}, -- configuration space upper right
   cells = {128}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   cflFrac = 1.0, -- CFL "fraction". Usually 1.0
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {8}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions
   writeGhost = true,
   -- electrons
   elc = Plasma.VlasovSpecies {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-vCells/2*dv},
      upper = {vCells/2*dv},
      cells = {vCells},
      decompCuts = {1},
      externalBC = "wall_16",
      -- initial conditions
      init = function (t, xn)
         local x, v = xn[1], xn[2]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
	    idx = 1280
	 else
	    idx = math.abs(x) * 10 / lambda_D
	    idx = math.floor(idx + 0.5) + 1
	 end
	 local n, u = n_e0[idx]*n_e, u_0[idx]*uB*sign
         return maxwellian(n, u, vth_e, v)
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.VlasovSpecies.bcReflect,
              Plasma.VlasovSpecies.bcExternal },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },

   ion = Plasma.VlasovSpecies {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-6.0*uB},
      upper = {6.0*uB},
      cells = {vCells},
      decompCuts = {1},
      -- initial conditions
      init = function (t, xn)
         local x, v = xn[1], xn[2]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
	    idx = 1280
	 else
	    idx = math.abs(x) * 10 / lambda_D
	    idx = math.floor(idx + 0.5) + 1
	 end
	 local n, u = n_i0[idx]*n_i, u_0[idx]*uB*sign
         return maxwellian(n, u, vth_i, v)
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.VlasovSpecies.bcReflect,
              Plasma.VlasovSpecies.bcAbsorb },
      diagnosticMoments = { "M0", "M1i", "M2" },
   },
   
   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = epsilon_0, mu0 = mu_0,
      init = function (t, xn)
         local x =  xn[1]
	 local idx = 1
	 local sign = math.abs(x)/x
	 if math.abs(x) > 128.0*lambda_D then
	    idx = 1280
	 else
	    idx = math.abs(x) * 10 / lambda_D
	    idx = math.floor(idx + 0.5) + 1
	 end
	 local E = E_0[idx]*elemCharge*lambda_D/T_e*sign
         return E, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.MaxwellField.bcSymmetry,
              Plasma.MaxwellField.bcReflect },
   },
}
-- run application
sim:run()
