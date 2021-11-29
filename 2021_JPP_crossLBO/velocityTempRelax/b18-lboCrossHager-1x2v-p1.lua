-- Gkyl ------------------------------------------------------------------------
--
--
local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"

-- This test initializes Maxwellian electrons and ions with different
-- bulk velocity and temperature and collides them.
-- Intended to reproduce figure 3 of R. Hager et al., JCP 315 (2016) 644–660.

-- Universal parameters.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE  -- (C).
hbar     = Constants.PLANCKS_CONSTANT_H/(2.*math.pi)        -- (J s).
me, mp   = Constants.ELECTRON_MASS, Constants.PROTON_MASS   -- (kg).

B0  = 1.0    -- Magnetic field amplitude [T].
n0  = 7e19   -- Number density [1/m^3].

mi     = 2.014*mp   -- Deuterium ions.
qe, qi = -eV, eV    -- Charge.
ne, ni = n0, n0     -- Density.
Te0    = 300*eV     -- Electron temperature.
Ti0    = 200*eV     -- Ion temperature.
alpha  = 1.3        -- Ratio of perpendicular to parallel temperature.

-- Parallel and perpendicular temperatures.
TePar0, TePerp0 = Te0, alpha*Te0
TiPar0, TiPerp0 = Ti0, alpha*Ti0

vti, vte = math.sqrt(Ti0/mi), math.sqrt(Te0/me)   -- Thermal speeds.

---- Bulk flow speed along B-field in terms of reference temperatures.
---- Plots indicate the possibility of these flows being computed in
---- of the true temperature of the IC, so we compute them that way later.
--uPare = 0.5*math.sqrt(me/mi)*vte
--uPari = 50.0*(me/mi)*vti

local function coulombLog_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   -- Coulomb logarithm in Gkeyll (see online documentation):
   local m_sr, u_sr = ms*mr/(ms+mr), math.sqrt(3.*vts^2+3.*vtr^2)
   local omega_ps, omega_pr = math.sqrt(ns*qs^2/(ms*eps0)), math.sqrt(nr*qr^2/(mr*eps0))
   local omega_cs, omega_cr = math.abs(qs*B0/ms), math.abs(qr*B0/mr)
   local rMax = 1./math.sqrt((omega_ps^2+omega_cs^2)/(vts^2+3.*vts^2)+(omega_pr^2+omega_cr^2)/(vtr^2+3.*vts^2))
   local rMin = math.max(math.abs(qs*qr)/(4.*math.pi*eps0*m_sr*u_sr^2), hbar/(2.*math.exp(0.5)*u_sr*m_sr))
   return 0.5*math.log(1. + (rMax/rMin)^2)
end
local function nu_ss(qs, ms, ns, us, vts) -- Like-species collision frequency.
   local logLambda = coulombLog_sr(qs, qs, ms, ms, ns, ns, vts, vts)
   return (1./math.sqrt(2))*(qs^4)*ns*logLambda/(3.*((2.*math.pi)^(3./2.))*(eps0^2)*(ms^2)*((vts^2)^(3./2.)))
end
local function nuM_sr(qs, qr, ms, mr, ns, nr, us, ur, vts, vtr)
   -- LBO-EM collision frequency.
   local logLambda = 0.5*( coulombLog_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
                          +coulombLog_sr(qr, qs, mr, ms, nr, ns, vtr, vts) )
   return 2.*(ms+mr)*((qs*qr)^2)*nr*logLambda/(3.*((2.*math.pi)^(3./2.))*(eps0^2)*ms^2*mr*((vts^2+vtr^2)^(3./2.)))
end
-- The ICs give a slightly different T than the reported reference temperature.
-- If Gkeyll were computing nu(x,t) from scratch that would not be a problem, but
-- since we are going to use normNu we need to use the true temperature of the ICs
-- in order to compute normNu:
TeIC, TiIC   = (2.*TePerp0+TePar0)/3., (2.*TiPerp0+TiPar0)/3. 
vtiIC, vteIC = math.sqrt(TiIC/mi), math.sqrt(TeIC/me)
-- Bulk flow speed along B-field.
uPare = 0.5*math.sqrt(me/mi)*vteIC
uPari = 50.0*(me/mi)*vtiIC
nu_ii = nu_ss(qi, mi, ni, uPari, vtiIC)
nu_ee = nu_ss(qe, me, ne, uPare, vteIC)
nu_ie = nuM_sr(qi, qe, mi, me, ni, ne, uPari, uPare, vtiIC, vteIC)
nu_ei = nuM_sr(qe, qi, me, mi, ne, ni, uPare, uPari, vteIC, vtiIC)
print(' Collision frequencies: ')
print('          nu_ii: ', nu_ii)
print('          nu_ee: ', nu_ee)
print('   LBO-EM nu_ie: ', nu_ie)
print('   LBO-EM nu_ei: ', nu_ei)
print(' ')
-- Normalize collision frequencies:
normNu_ee, normNu_ei = nu_ee*((2.*vteIC^2)^(3/2))/ne, nu_ei*((vteIC^2+vtiIC^2)^(3/2))/ni
normNu_ii, normNu_ie = nu_ii*((2.*vtiIC^2)^(3/2))/ni, nu_ie*((vtiIC^2+vteIC^2)^(3/2))/ne

-- bi-Maxwellian distribution with drift u and temperature T.
local function biMaxwellian(n, upar, T, m, B, alpha, x, vpar, mu)
   local vExp = -(0.5*m*((vpar-upar)^2+(2.*mu*B/m)/alpha)/T)
   return (n/(alpha*(2.0*math.pi*(T/m))^(3./2.)))*math.exp(vExp)
end

plasmaApp = Plasma.App {
   logToFile = true,
   
   tEnd        = tEndFrac/nu_ii,      -- End time.
   nFrame      = numFrames,       -- Number of frames to write.
   lower       = {-2.},         -- Configuration space lower coordinate.
   upper       = { 2.},         -- Configuration space upper coordinate.
   cells       = {1},            -- Configuration space cells.
   basis       = "serendipity",   -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,               -- Polynomial order.
   timeStepper = "rk3",           -- One of "rk2", "rk3" or "rk3s4".
   cflFrac     = 1.0,
   calcIntQuantEvery = 1./(5.*numFrames),
   groupDiagnostics = true,
   
   -- Decomposition for configuration space.
   decompCuts = {1},              -- Cuts in each configuration direction.
   useShared  = false,            -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1},            -- Periodic directions.

   -- Neutral species with a rectangular/square IC.
   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {-5.*vte, 0.0},
      upper = { 5.*vte, me*(5.*vte)^2/(2.*B0)},
      cells = {16, 16},
      -- Initial conditions.
      init = Plasma.FunctionProjection {
         func = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return biMaxwellian(n0, uPare, Te0, me, B0, alpha, x, vpar, mu)
         end, 
      },
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnostics = { "M0", "Upar", "Tpar", "Tperp", "intM0", "intM1", "intM2"},
      nDistFuncFrame = 1,
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {"elc", "ion" },
	 normNu = {normNu_ee, normNu_ei},
         lboType = "LBO-EM",
      },
   },

   -- Neutral species with a bump in the tail.
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {-5.*vti, 0.0},
      upper = { 5.*vti, mi*(5.*vti)^2/(2.*B0)},
      cells = {16, 16},
      -- Initial conditions.
      init = Plasma.FunctionProjection {
         func = function (t, xn)
            local x, vpar, mu = xn[1], xn[2], xn[3]
            return biMaxwellian(n0, uPari, Ti0, mi, B0, alpha, x, vpar, mu)
         end, 
      },
      -- Evolve species?
      evolve = true,
      evolveCollisionless = false,
      -- Diagnostic moments.
      diagnostics = { "M0", "Upar", "Tpar", "Tperp", "intM0", "intM1", "intM2"},
      nDistFuncFrame = 1,
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith = {"ion", "elc"},
	 normNu = {normNu_ii, normNu_ie},
         lboType = "LBO-EM",
      },
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = false,    -- Evolve fields?
      kperp2 = 0.0 
   },
   
   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      bmag = function (t, xn) return B0 end,
      evolve = false,  -- Geometry is not time-dependent.
   },

}
-- Run application.
plasmaApp:run()
