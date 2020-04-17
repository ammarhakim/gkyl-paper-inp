-- Gkyl ------------------------------------------------------------------------
--
-- Helimak simulation.
--
---------------------------------------------------------------------------

local Plasma    = require("App.PlasmaOnCartGrid").Gyrokinetic
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"

-- Physical parameters.
eps0 = Constants.EPSILON0             -- Permittivity of vacuum (F/m).
eV   = Constants.ELEMENTARY_CHARGE
qe   = -eV                            -- Electron charge (C).
qi   = eV                             -- Ion charge (C).
mi   = 39.948*Constants.PROTON_MASS   -- Argon ion mass (kg).
me   = mi/400                         -- Electron mass (kg).

Te0 = 10*eV   -- Electron temperature (J).
Ti0 = 1*eV    -- Ion temperature (J).
R0  = 1.1     -- Major radius (m).
B0  = 0.1     -- Magnetic field at R=R0 (T).
n0  = 1e16    -- Reference density (1/m^3).
np  = 9*n0

vti      = math.sqrt(Ti0/mi)    -- Ion thermal speed (m/s).
vte      = math.sqrt(Te0/me)    -- Electron thermal speed (m/s).
c_s      = math.sqrt(Te0/mi)    -- Sound speed (m/s).
omega_ci = math.abs(qi*B0/mi)   -- Ion cyclotron frequency (rad/s).
rho_s    = c_s/omega_ci         -- Ion sound gyroradius (m).

-- Simulation box size (m).
Lx = 50*rho_s
Ly = 17*rho_s
Lz = 40

-- Source parameters.
Ts             = 5.0/3.0*(Te0+Ti0)*eV
c_ss           = math.sqrt(Ts/mi)
S0             = 3e19
x_nSource      = 1.0    -- Source start coordinate (m).
lambda_nSource = 0.01   -- Characteristic length scale of density and temperature (m).

-- Collision parameters.
nuFrac       = 1.0
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc        = nuFrac*logLambdaElc*eV^4*n0
              /(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))

logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon        = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))

-- Initial profiles.
initDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return math.exp(-(x-x_nSource)^2/(2*lambda_nSource^2))+0.1
end
initTemperature = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 1
end
-- Source profiles.
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 4*S0*math.exp(-(x-x_nSource)^2/(2*lambda_nSource^2))+0.001
end
sourceTemperatureElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 12.0/3.0*Te0
end
sourceTemperatureIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 5.0/3.0*Ti0
end

-- Initialize a random seed for initial conditions.
-- Will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd              = 16e-3,                     -- End time.
   nFrame            = 1600,                      -- Number of output frames.
   lower             = {R0 - Lx/2, -Ly/2, -Lz/2}, -- Configuration space lower left.
   upper             = {R0 + Lx/2, Ly/2, Lz/2},   -- Configuration space upper right.
   cells             = {48, 24, 16},              -- Configuration space cells.
   basis             = "serendipity",             -- One of "serendipity" or "maximal-order".
   polyOrder         = 1,                         -- Polynomial order.
   timeStepper       = "rk3",                     -- One of "rk2" or "rk3".
   cflFrac           = 0.2,
   restartFrameEvery = 0.01,

   -- Decomposition for configuration space.
   decompCuts = {8, 8, 4}, -- Cuts in each configuration direction.
   useShared  = false,     -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {2},   -- Periodic in y only.

   -- Gyrokinetic electrons.
   electron = Plasma.Species {
      charge = qe,
      mass   = me,
      -- Velocity space grid.
      lower      = {-4*vte, 0},
      upper      = { 4*vte, 12*me*vte^2/(4*B0)},
      cells      = {10, 5},
      -- Initial conditions.
      init = {"maxwellian",
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local effectiveSource = np*initDensity(t,{x,y,0})
                 local perturb = 1e-3*(math.random()-0.5)*2
                 if (1 - z^2/(Lz/2)^2) < 0 then
                    return 0.5*effectiveSource*(1+perturb)
                 else
		    return effectiveSource*(1+math.sqrt(1-(2*z/Lz)^2))/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 return Te0*initTemperature(t,{x,y,z})
              end,
             },
      coll = Plasma.LBOCollisions {
	  frequencies = {nuElc,},
	  collideWith = {'electron'},
      },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperatureElc},
      evolve = true, -- Evolve species?
      diagnosticIntegratedMoments = { "intM0", "intM2"},
      diagnosticMoments           = {"GkM0", "GkM1"},
      randomseed = randomseed,
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath,   Plasma.Species.bcSheath},
   },

   -- Gyrokinetic ions.
   ion = Plasma.Species {
      charge = qi,
      mass   = mi,
      -- Velocity space grid.
      lower = {-6*vti, 0},
      upper = {6*vti, 12*mi*vti^2/(2*B0)},
      cells = {10, 5},
      -- Initial conditions.
      init = {"maxwellian",
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local effectiveSource   = np*initDensity(t,{x,y,0})
                 local perturb           = 1e-3*(math.random()-0.5)*2
                 if (1 - z^2/(Lz/2)^2) < 0 then
                    return 0.5*effectiveSource*(1+perturb)
		 else
		    return effectiveSource*(1+math.sqrt(1-(2*z/Lz)^2))/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
		 return Ti0*initTemperature(t,{x,y,z})
              end,
              driftSpeed = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 if z == 0 then
                    return 0
                 else
                    return math.sqrt(3.0)/4.0*math.sqrt(Ts*eV/mi)*(1 - math.sqrt(1 - (2*z/Lz)^2))*Lz/z 
                 end
              end,
             },
      coll = Plasma.LBOCollisions {
         frequencies = {nuIon, },
	 collideWith = {'ion'},
      },
      source = {"maxwellian", density = sourceDensity, temperature = sourceTemperatureIon},
      evolve = true, -- Evolve species?
      diagnosticMoments           = {"GkM0", "GkM1"},
      diagnosticIntegratedMoments = { "intM0", "intM2"},
      randomseed = randomseed,
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath,   Plasma.Species.bcSheath},
   },

   -- Field solver.
   field = Plasma.Field {
      -- Dirichlet in x.
      phiBcLeft  = { T ="D", V = 0.0},
      phiBcRight = { T ="D", V = 0.0},
      -- Periodic in y --
      -- No bc in z.
      phiBcBack  = { T ="N", V = 0.0},
      phiBcFront = { T ="N", V = 0.0},
      evolve     = true,   -- Evolve fields?
   },

   -- Magnetic geometry .
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R0/x
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
