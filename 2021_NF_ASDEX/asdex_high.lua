-- Plasma ------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"

-- scan parameters
Tfac = 1
nfac = 1*math.sqrt(Tfac)

-- Universal constant parameters.
eps0         = Constants.EPSILON0
eV           = Constants.ELEMENTARY_CHARGE
qe           = -eV
qi           = eV
me           = Constants.ELECTRON_MASS

-- Plasma parameters.
mi           = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0          = 85*eV*Tfac
Ti0          = 85*eV*Tfac 
n0           = 2e19*nfac  -- [1/m^3]

-- Geometry and magnetic field.
B_axis       = 2.45   -- [T]
R0           = 1.65  -- [m]
a0           = 0.5   -- [m]
R            = R0 + a0
B0           = B_axis*(R0/R) -- [T]
Lpol         = 2.4 -- [m]

-- Parameters for collisions.
nuFrac = 1.0
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))
-- Electron-Ion collision frequency.
nuElcIon     = nuElc/1.96
nuIonElc     = me*nuElcIon/mi -- Ion-electron collision frequency.

-- Derived parameters
vti      = math.sqrt(Ti0/mi)
vte  	 = math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci

-- Box size.
--Lx = 50*rho_s
Ly = 100*rho_s
Lz = 10 -- [m]

-- Source parameters.
P_SOL        = 0.8e6*nfac*Tfac -- [W] 
P_src        = P_SOL*Ly*Lz/(2*math.pi*R*Lpol)
xSource      = R - 0.003 -- [m], source start coordinate
lambdaSource = 0.001 -- [m], characteristic length scale of density and temperature

-- Source profiles.
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local sourceFloor = 1e-10
   if x < xSource + 3*lambdaSource then
      sourceFloor = 1e-2
   end
   if math.abs(z) < Lz/4 then
      return math.max(math.exp(-(x-xSource)^2/(2*lambdaSource)^2), sourceFloor)
   else
      return 1e-40
   end
end
sourceTemperature = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < xSource + 3*lambdaSource then
      return 90*eV*Tfac
   else
      return 5*eV*Tfac
   end
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 4e-3,                     -- End time.
   nFrame      = 4000,                     -- Number of output frames.
   --lower       = {R - Lx/2-.004, -Ly/2, -Lz/2}, -- Configuration space lower left.
   lower       = {2.142, -Ly/2, -Lz/2},
   --upper       = {R + Lx/2, Ly/2, Lz/2},   -- Configuration space upper right.
   upper       = {2.187, Ly/2, Lz/2},
   cells       = {16, 32, 10},              -- Configuration space cells.
   basis       = "serendipity",            -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                        -- Polynomial order.
   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.9,
   restartFrameEvery = 0.001,

   -- Decomposition for configuration space.
   decompCuts = {8, 16, 2}, -- Cuts in each configuration direction.
   useShared = false,      -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {2},     -- Periodic in y only.

   -- Gyrokinetic electrons.
   electron = Plasma.Species {
      charge = qe,
      mass  = me,
      lower = {-4*vte, 0},
      upper = {4*vte, 12*me*vte^2/(2*B0)},
      cells = {10, 5},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local Ls              = Lz/4
                 local floor = 0.01
                 local effectiveSource = math.max(sourceDensity(t,{x,y,0}), floor)
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb         = 1e-3*(math.random()-0.5)*2.0
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if (x < xSource + 3*lambdaSource) then 
                    return 85*eV*Tfac
                 else 
                    return 5*eV*Tfac
                 end
              end,
              scaleWithSourcePower = true,
      },
      coll   = Plasma.LBOCollisions { 
         collideWith = {'electron', 'ion'},
         frequencies = {nuElc, nuElcIon},
         --nuFrac = nuFrac,
      },
      source = Plasma.MaxwellianProjection {
                density = sourceDensity,
                temperature = sourceTemperature,
                power = P_src/2,
                isSource = true,
      },
      evolve = true, -- Evolve species?
      diagnosticMoments = {"GkM0", "GkUpar", "GkTemp", "GkBeta"}, 
      diagnosticIntegratedMoments = {"intM0", "intKE", "intHE", "intL2"},
      diagnosticBoundaryFluxMoments = {"GkM0", "GkM1", "GkHamilEnergy"},
      diagnosticIntegratedBoundaryFluxMoments = {"intM0", "intM1", "intHE"},
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      charge = qi,
      mass   = mi,
      -- Velocity space grid.
      lower = {-4*vti, 0},
      upper = {4*vti, 12*mi*vti^2/(2*B0)},
      cells = {10, 5},
      decompCuts = {1, 1},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x, y, z         = xn[1], xn[2], xn[3]
                 local Ls              = Lz/4
                 local floor = 0.01
                 local effectiveSource = math.max(sourceDensity(t,{x,y,0}), floor)
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi)
                 local nPeak           = 4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
                 local perturb         = 1e-3*(math.random()-0.5)*2.0
                 if math.abs(z) <= Ls then
                    return nPeak*(1+math.sqrt(1-(z/Ls)^2))/2*(1+perturb)
                 else
                    return nPeak/2*(1+perturb)
                 end
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 if x < xSource + 3*lambdaSource then 
                    return 85*eV*Tfac
                 else 
                    return 5*eV*Tfac
                 end
              end,
              driftSpeed = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local Te
                 if x < xSource + 3*lambdaSource then 
                    Te = 85*eV
                 else 
                    Te = 5*eV
                 end
                 if math.abs(z) <= Lz/4 then
		    return z/(Lz/4)*math.sqrt(Te/mi)
                 else
		    return z/math.abs(z)*math.sqrt(Te/mi)
                 end
              end,
              scaleWithSourcePower = true,
      },
      coll   = Plasma.LBOCollisions { 
         collideWith = {'ion', 'electron'},
         frequencies = {nuIon, nuIonElc},
         --nuFrac = nuFrac,
      },
      source = Plasma.MaxwellianProjection {
                density = sourceDensity,
                temperature = sourceTemperature,
                power = P_src/2,
                isSource = true,
      },
      evolve = true, -- Evolve species?
      diagnosticMoments = {"GkM0", "GkUpar", "GkTemp", "GkBeta"}, 
      diagnosticIntegratedMoments = {"intM0", "intKE", "intHE", "intL2"},
      diagnosticBoundaryFluxMoments = {"GkM0", "GkM1", "GkHamilEnergy"},
      diagnosticIntegratedBoundaryFluxMoments = {"intM0", "intM1", "intHE"},
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},
   },

   -- Field solver.
   field = Plasma.Field {
      -- Dirichlet in x.
      phiBcLeft  = { T ="D", V = 0.0},
      phiBcRight = { T ="D", V = 0.0},
      aparBcLeft  = { T ="D", V = 0.0},
      aparBcRight = { T ="D", V = 0.0},
      -- Periodic in y. --
      -- No bc in z.
      phiBcBack  = { T ="N", V = 0.0},
      phiBcFront = { T ="N", V = 0.0},
      evolve     = true, -- Evolve fields?
      isElectromagnetic = true,
   },

   -- Magnetic geometry.
   funcField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R/x
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- run application
plasmaApp:run()
