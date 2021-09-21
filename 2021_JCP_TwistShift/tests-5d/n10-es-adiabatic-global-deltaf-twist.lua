-- Input file based on Gorler et al, Phys Plasmas 23 (2016)
-- ITG with adiabatic electron case, see section IV.A of paper
-- 'Local' geometry (narrow domain in x)

-- Plasma ------------------------------------------------------------------------
local Plasma = require ("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local math = require("sci.math").generic

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qi = eV
qe = -eV

-- _reference parameters (see Table II of Gorler et al)
m_ref = 2.0*Constants.PROTON_MASS -- (deuterium ions)
n_ref = 4.66e19 -- m^-3
T_ref = 2.14e3*eV 
B_ref = 2.0 -- T
L_ref = 1.67 -- m
beta_ref = 0.0101
rhostar = 1/180.2  -- = rho_s/a

-- parameters (see Table I of Gorler et al)
a = 0.36*L_ref
r0 = 0.5*a
R0 = L_ref
T0 = T_ref
kap_T = 6.96
width_T = 0.3
n0 = n_ref
kap_n = 2.23
width_n = 0.3
mi = m_ref
me = 5.44617e-4*m_ref
toroidal_mode_number = 10  -- = n_0 (scan parameter)

-- profiles
function qprofile(r)
   return 2.52*(r/a)^2 - 0.16*(r/a) + 0.86
end

function nTprofile(r, A0, kap, width)
   return A0*math.exp(-kap*width*a/L_ref*math.tanh((r-r0)/(width*a)))
end

-- derived parameters
c_s = math.sqrt(T_ref/m_ref)
vti = c_s
vte = math.sqrt(T_ref/me)
omega_ref = qi*B_ref/(m_ref)
rho_s = c_s/omega_ref
q0 = qprofile(r0)
ky = toroidal_mode_number*q0/r0

-- grid parameters
Ly = 2*math.pi/ky
Lx = 0.8*a
nperiod = 1 -- number of 2*pi periods along field line
Lz = 2*math.pi*nperiod
N_VPAR, N_MU = 16, 8
VPAR_UPPER = 3*vti
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = 9*T_ref/B_ref

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 20*R0/vti, -- end time
   nFrame = 100, -- number of output frames
   lower = {r0-Lx/2, -Ly/2, -Lz/2}, -- configuration space lower left
   upper = {r0+Lx/2,  Ly/2, Lz/2}, -- configuration space upper right
   cells = {96, 16, 16*nperiod}, -- configuration space cells
   mapc2p = function(xc)
      local x, y, z = xc[1], xc[2], xc[3]
      local q = qprofile(x)

      local phi = q*z - q0/r0*y
      -- local theta = z 
      local eps = x/R0
      local costheta = (math.cos(z) - eps)/(1 - eps*math.cos(z))
      local sintheta = math.sqrt(1-eps^2)/(1 - eps*math.cos(z))*math.sin(z)

      -- map to cylindrical (R, Z, phi) coordinates
      local R = R0 + x*costheta
      local Z = x*sintheta

      -- map to Cartesian (X, Y, Z) coordinates
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)

      return X, Y, Z
   end,
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 0.5,

   -- decomposition for configuration space
   decompCuts = {24, 1, 8}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {2}, -- periodic directions
   --deltaF = true,
   groupDiagnostics = true,

   -- gyrokinetic ions
   ion = Plasma.Species {
      charge = qi,
      mass = mi,
      -- velocity space grid
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      decompCuts = {1, 1},
      -- initial conditions
      background = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x = xn[1]
                 return nTprofile(x, n_ref, kap_n, width_n)
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 return nTprofile(x, T_ref, kap_T, width_T)
              end,
              exactScaleM012 = true,
      },
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local perturb = 1e-10*rho_s/a*math.cos(toroidal_mode_number*qprofile(r0)/r0*y) -- should this be q(r0) or q(r)?

                 -- filter perturbation at domain ends in z so that fluctuations -> 0
                 local zfilter = 1.0
                 if math.abs(z-nperiod*math.pi)<math.pi then
                    local w = (z - nperiod*math.pi)/(math.pi)
                    zfilter = 2*w^2/(1+w^4)
                 end
                 if math.abs(z+nperiod*math.pi)<math.pi then
                    local w = (z + nperiod*math.pi)/(math.pi)
                    zfilter = 2*w^2/(1+w^4)
                 end

                 return 0*nTprofile(x, n_ref, kap_n, width_n) + n_ref*perturb*zfilter
              end,
              driftSpeed = 0.0,
              temperature = function (t, xn)
                 local x = xn[1]
                 return nTprofile(x, T_ref, kap_T, width_T)
              end,
              --exactScaleM012 = true,
      },
      evolve = true, 
      diagnostics = {"M0", "M1", "M2"}, --, perturbed=true}, 
      bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}},
      bcz = {Plasma.TwistShiftBC{shiftFunction=function(t,xn) return r0/q0*qprofile(xn[1])*Lz end},
             Plasma.TwistShiftBC{shiftFunction=function(t,xn) return r0/q0*qprofile(xn[1])*Lz end}},
      deltafGK = true,
      deltafLinear = true,
   },

   -- adiabatic electrons
   adiabaticElectron = Plasma.AdiabaticSpecies {
      charge = qe,
      mass = me,
      temp = T_ref,
      -- initial conditions
      init = function (t, xn)
         local x = xn[1]
         return 0*nTprofile(x, n_ref, kap_n, width_n)
      end,
      evolve = false, 
   },

   -- field solver
   field = Plasma.Field {
      evolve = true, -- evolve fields?
      bcLowerPhi  = {{T = "N", V = 0.0}, {T = "P"}, {T = "N", V = 0.0}},
      bcUpperPhi  = {{T = "N", V = 0.0}, {T = "P"}, {T = "N", V = 0.0}},
   },

   -- magnetic geometry 
   extField = Plasma.Geometry {
      bmag = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]	
         local q = qprofile(x)
         -- local theta = z 
         local eps = x/R0
         local costheta = (math.cos(z) - eps)/(1 - eps*math.cos(z))

         -- map to cylindrical (R, Z, phi) coordinates
         local R = R0 + x*costheta

         local Bt = B_ref*R0/R
         local Bp = Bt*eps/q/math.sqrt(1-eps^2)

         return math.sqrt(Bt^2 + Bp^2)
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
