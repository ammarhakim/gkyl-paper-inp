-- Gkyl ------------------------------------------------------------------------
--
-- Collide two 1x1v p=1 Vlasov species with different mass, mean flow and/or
-- temperature with the multispecies LBO-G operator.
--
local Plasma    = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Constants = require "Lib.Constants"

-- Universal parameters (SI units).
local eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
local hbar     = Constants.PLANCKS_CONSTANT_H/(2.*math.pi) -- (J s).
local qe, qi   = -eV, eV 
local mp = Constants.PROTON_MASS
local me = Constants.ELECTRON_MASS

local n0  = 7e19          -- Number density (1/m^3).
local Te0 = 40*eV         -- Electron temperature (J).
local Ti0 = 80*eV         -- Ion temperature (J).
local mi  = mp            -- Ion mass (kg).
local massRatio = mi/me   -- Ion to electron mass ratio.
local B0  = 0.            -- External magnetic field amplitude (T).

local me = mi/massRatio   -- Electron mass (kg).
local vte, vti = math.sqrt(Te0/me), math.sqrt(Ti0/mi)   -- Thermal speeds (m/s).

local polyOrder    = 1
local vMini, vMaxi = -5.*vti, 5.*vti   -- Minimum and maximum ion velocity. 
local vMine, vMaxe = -5.*vte, 5.*vte   -- Minimum and maximum electron velocity. 
local Nv           = 8                 -- Number of cells in one dim of velocity space.

-- Parameters defining the initial condition.
local A, delta       = n0, 0.5
local kve, kvi       = (2.*math.pi/(2.*vte)),  (2.*math.pi/(2.*vti))
local betae, betai   = {0.5*vte, -0.5*vte}, {1.5*vti, 1.5*vti}
local sigmae, sigmai = vte, vti

-- Initial distribution function.
local function initDistF(A, delta, kv, beta, sigma, v)
   local vMag = math.sqrt(v[1]^2+v[2]^2)
   return A*(1.+delta*math.cos(kv*vMag))*math.exp(-((v[1]-beta[1])^2+(v[2]-beta[2])^2)/(2.*(sigma^2)))
end

local function logLambda_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   -- Coulomb logarithm in Gkeyll (see online documentation):
   local m_sr, u_sr = ms*mr/(ms+mr), math.sqrt(3.*vts^2+3.*vtr^2)
   local omega_ps, omega_pr = math.sqrt(ns*qs^2/(ms*eps0)), math.sqrt(nr*qr^2/(mr*eps0))
   local omega_cs, omega_cr = math.abs(qs*B0/ms), math.abs(qr*B0/mr)
   local rMax = 1./math.sqrt((omega_ps^2+omega_cs^2)/(vts^2+3.*vts^2)+(omega_pr^2+omega_cr^2)/(vtr^2+3.*vts^2))
   local rMin = math.max(math.abs(qs*qr)/(4.*math.pi*eps0*m_sr*u_sr^2), hbar/(2.*math.exp(0.5)*u_sr*m_sr))
   return 0.5*math.log(1. + (rMax/rMin)^2)
end
local function nuM_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   -- LBO-EM collision frequency.
   local coulombLog = logLambda_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   return 2.*(ms+mr)*((qs*qr)^2)*nr*coulombLog/(3.*((2.*math.pi)^(3./2.))*(eps0^2)*ms^2*mr*((vts^2+vtr^2)^(3./2.)))
end
local function nuT_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   -- LBO-ET collision frequency.
   local coulombLog = logLambda_sr(qs, qr, ms, mr, ns, nr, vts, vtr)
   return 2.*((qs*qr)^2)*nr*coulombLog/(3.*((2.*math.pi)^(3./2.))*(eps0^2)*ms*mr*((vts^2+vtr^2)^(3./2.)))
end

local nuElcIon = nuM_sr(qe, qi, me, mi, n0, n0, vte, vti)
local nuIonElc = nuM_sr(qi, qe, mi, me, n0, n0, vti, vte)
local nuElcElc = nuM_sr(qe, qe, me, me, n0, n0, vte, vte)
local nuIonIon = nuM_sr(qi, qi, mi, mi, n0, n0, vti, vti)

print(" Collision frequencies (1/s):")
print(string.format("   nu_ei = %g", nuElcIon))
print(string.format("   nu_ie = %g", nuIonElc))
print(string.format("   nu_ee = %g", nuElcElc))
print(string.format("   nu_ii = %g", nuIonIon))

-- Estimate the time step and set the simulation end at 1e4 time steps.
local dv         = (vMaxe-vMine)/Nv
local vDim       = #betae
local lambdaMax  = 4.*nuElcIon*((polyOrder+1)^2)*(vte^2)*(1./(dv^2))*vDim
local dtApprox   = 1./lambdaMax
print(string.format(" Approximate time step (s): %g", dtApprox))

vlasovApp = Plasma.App {
   tEnd        = 1.e4*dtApprox,    -- End time.
   nFrame      = 10,               -- Number of frames to write.
   lower       = {-1.0},           -- Configuration space lower left.
   upper       = { 1.0},           -- Configuration space upper right.
   cells       = {1} ,             -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = polyOrder,        -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2", "rk3" or "rk3s4".
   maximumDt   = dtApprox,
   cflFrac     = 1e6,              -- To force dt=maximumDt.

   decompCuts = {1},          -- Cuts in each configuration direction.
   useShared  = true,        -- If to use shared memory.
   groupDiagnostics = true,   -- Group some diagnostics into a single file.

   -- Boundary conditions for configuration space.
   periodicDirs = {1}, -- Periodic directions.

   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {vMini, vMini},
      upper = {vMaxi, vMaxi},
      cells = {Nv, Nv},
      init = function(t, xn)    -- Initial conditions.
	 x, v = xn[1], {xn[2], xn[3]}
         return initDistF(A, delta, kvi, betai, sigmai, v)
      end,
      evolve = true,
      evolveCollisionless = false,  -- No collisionless terms.
      diagnosticIntegratedMoments = {"intM0", "intM1i", "intM2"},
      nDistFuncFrame = 1,  -- Output a single distribution function.
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith  = { "elc" },
         frequencies  = { nuIonElc },
         -- Optional arguments:
         lboType      = "LBO-EM",
      },
   },

   elc = Plasma.Species {
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {vMine, vMine},
      upper = {vMaxe, vMaxe},
      cells = {Nv, Nv},
      init = function(t, xn)    -- Initial conditions.
	 x, v = xn[1], {xn[2], xn[3]}
         return initDistF(A, delta, kve, betae, sigmae, v)
      end,
      evolve = true,
      evolveCollisionless = false,  -- No collisionless terms.
      diagnosticIntegratedMoments = {"intM0", "intM1i", "intM2"},
      nDistFuncFrame = 1,  -- Output a single distribution function.
      -- Collisions.
      coll = Plasma.LBOCollisions {
         collideWith  = { "ion"},
         frequencies  = { nuElcIon },
         -- Optional arguments:
         lboType      = "LBO-EM",
      },
   },

}
-- Run application.
vlasovApp:run()
