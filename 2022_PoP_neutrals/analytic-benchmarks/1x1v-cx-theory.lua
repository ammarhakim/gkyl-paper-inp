-- Gkyl ----------------------------------------------------------------
-- Basic sheath simulation ---------------------------------------------
-- To test neutral model in 1x1v with static plasma and charge exchange.
-- Results produced with Gkeyll changeset [d9eafc9f1234].
------------------------------------------------------------------------

local Plasma = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"

function sech(x)
   return 2*math.exp(x)/(math.exp(2*x)+1)
end


nfac = 1.0
-- Universal parameters.
eps0 = Constants.EPSILON0
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV
qi   =  eV
mi   = Constants.PROTON_MASS   -- Hydrogen ions.
me   = Constants.ELECTRON_MASS -- True electron mass.
B0  = 0.5                      -- Magnetic field amplitude [T].

n0  = nfac*5e18                -- Number density [1/m^3].

Te0 = 20*eV                    -- Electron temperature.
Ti0 = 20*eV                    -- Ion temperature.

-- Derived physical parameters
-- set kmin*rho_s = 0.2, used for Field table
cs = math.sqrt(Te0/mi)
omega_ci = eV*B0/mi
rho_s = cs/omega_ci
kmin = 0.2/rho_s

-- Thermal speeds.
vti = math.sqrt(Ti0/mi)
vte = math.sqrt(Te0/me)

-- Connection length
Lx = 40 --[m]

-- Neutral source parameters
sourceDensityNeut = function (t,xn)
   local x = xn[1]
   local x0 = 1
   local S0n = 8e-21*n0^2      -- to approximate recombination source
   return S0n
end

sim = Plasma.App {
   logToFile = false,

   tEnd        = 3*Lx/cs,          -- End time.
   nFrame      = 30,               -- Number of output frames.
   lower       = {-Lx/2},          -- Configuration space lower left.
   upper       = {Lx/2},           -- Configuration space upper right.
   cells       = {224},            -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   cflFrac     = 1,                -- CFL "fraction". Usually 1.0.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".

   -- Decomposition for configuration space.
   decompCuts = {4},   -- Cuts in each configuration direction.
   useShared  = false,  -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {}, -- Periodic directions.
      
    -- Electrons.
   elc = Plasma.Species {
      nDistFuncFrame = 10,
      charge = qe, mass = me,
      -- Velocity space grid.
      lower = {-4.0*vte},
      upper = {4.0*vte},
      cells = {16},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection{
         density = function (t, xn)
	    return n0
	 end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
	    if x <= 0 then
	       return -cs
	    else
	       return cs
	    end
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return Te0
         end,
      },
      evolve = false, -- Evolve species?
      -- Diagnostics
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
      				     "intM2" },
   },

   -- Ions
   ion = Plasma.Species {
      charge = qi, mass = mi,
      -- Velocity space grid.
      lower = {-4.0*vti},
      upper = {4.0*vti},
      cells = {16},
      decompCuts = {1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
	    return n0
	 end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2] 
	    if x <= 0 then
	       return -cs
	    else
	       return cs
	    end
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
	    return Ti0
         end,
      },
      evolve = false,
      -- Neutral interactions
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"neut"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      	 charge = qi,
	 vSigmaCX = 2.2e-14,
      },
      -- Diagnostics
      diagnosticMoments = { "GkM0", "GkM1", "GkM2", "GkUpar", "GkVtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1",
				     "intM2" },
   },

   neut = Plasma.Vlasov {
      charge = 0.0, mass = mi,
      -- Velocity space grid
      lower = {-4.0*vti},
      upper = {4.0*vti},
      cells = {32},
      decompCuts = {1},
      init = Plasma.VmMaxwellianProjection {
         density = function (t, xn)
            local x, vpar = xn[1], xn[2]
            local n_n = n0
	    local x0 = 0.2
	    local flr = 1e-6
	    if x <= 0 then
	       return n_n*(sech((-Lx/2-x)/x0)^2 + flr)
	    else	       
	       return n_n*(sech((Lx/2-x)/x0)^2 + flr)
	    end
	    -- return 0.1*n0
         end,
         driftSpeed = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return {0,0,0}
         end,
         temperature = function (t, xn)
            local x, vpar = xn[1], xn[2]
            return 2*eV
         end,
      },
      evolve = true,
      -- Boundary conditions
      bcx = {Plasma.Vlasov.bcRecycle, Plasma.Vlasov.bcRecycle},
      recycleTemp = 2*eV,
      recycleFrac = 1.0,
      recycleSpeed = cs,
      recycleFluxFac = 4.977492882411091e18,
      recycleIon = "ion",
      -- Neutral interactions
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"ion"},
      	 ions = "ion",
      	 neutrals = "neut",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "H",
      	 charge = 0,
	 vSigmaCX = 2.2e-14,
      },
      -- Diagnostics
      diagnosticMoments = { "M0", "u", "vtSq"},
      diagnosticIntegratedMoments = {"intM0", "intM1i",
      				     "intM2Flow", "intM2Thermal" },
      diagnosticBoundaryFluxMoments = {"M0"},
   },
   
   -- Field solver.
   field = Plasma.Field {
      bcLowerPhi = {{ T ="N", V = 0.0}},
      bcUpperPhi = {{ T ="N", V = 0.0}},
      evolve = true,
      kperp2 = kmin*kmin,
   },

      -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- Run application.
sim:run()
