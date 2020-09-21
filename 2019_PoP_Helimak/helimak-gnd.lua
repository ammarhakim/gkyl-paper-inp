-- Helimak simulation for Gkeyll v2.0, 
-- adapted from Lua input file for Gkeyll v1.0.
--
-- Plasma ------------------------------------------------------------------------
local Plasma = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi = require "Comm.Mpi"

-- physical parameters
eps0 = Constants.EPSILON0
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
mi = 39.948*Constants.PROTON_MASS -- (Argon ions)
me = mi/400                       -- reduced ion-to-electron mass ratio

Te0 = 10*eV 
Ti0 = 1*eV 
R0 = 1.1     -- [m]
B0 = 0.1     -- [T]
n0 = 1e16    -- [1/m^3]
np = 9*n0    -- [1/m^3] for initial density profile

-- derived parameters
vti = math.sqrt(Ti0/mi)        -- [m/s] ion thermal velocity
vte = math.sqrt(Te0/me)        -- [m/s] electron thermal velocity
c_s = math.sqrt(Te0/mi)        -- [m/s] ion sound speed
omega_ci = math.abs(qi*B0/mi)  -- [1/s] ion gyro-frequency 
rho_s = c_s/omega_ci           -- [m] ion gyro-radius

-- box size
Lx = 50*rho_s  -- [m]
Ly = 17*rho_s  -- [m]
Lz = 40        -- [m]

-- source parameters
Ts = 5.0/3.0*(Te0+Ti0)*eV
c_ss = math.sqrt(Ts/mi)
x_nSource = 1.0       -- [m], source start coordinate
lambda_nSource = 0.01 -- [m], characteristic length scale of density
P_src = 150           -- [W], source is scaled to input power

-- initial condition profiles
initDensity = function (t, xn)                                                            
   local x, y, z = xn[1], xn[2], xn[3]
   return math.exp(-(x-x_nSource)^2/(2*lambda_nSource^2))+0.1
end

-- source profiles
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return math.exp(-(x-x_nSource)^2/(2*lambda_nSource^2))+1e-3
end
sourceTemperatureElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 10.0/3.0*Te0
end
sourceTemperatureIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 5.0/3.0*Ti0
end

-- initialize a random seed for initial conditions
-- will be used for both ions and electrons
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 16e-3,                         -- end time
   nFrame = 1600,                        -- number of output frames
   lower = {R0 - Lx/2, -Ly/2, -Lz/2},    -- configuration space lower left
   upper = {R0 + Lx/2, Ly/2, Lz/2},      -- configuration space upper right
   cells = {48, 24, 16},                 -- configuration space cells
   basis = "serendipity",                -- one of "serendipity" or "maximal-order"
   polyOrder = 1,                        -- polynomial order
   timeStepper = "rk3",                  -- one of "rk2" or "rk3"
   cflFrac = 0.2,
   restartFrameEvery = 0.01,

   -- decomposition for configuration space
   decompCuts = {8, 8, 4},               -- cuts in each configuration direction
   useShared = false,                    -- shared memory is off

   -- boundary conditions for configuration space
   periodicDirs = {2}, -- periodic in y only

   -- gyrokinetic electrons
   electron = Plasma.Species {
      evolve = true, -- evolve species?
      charge = qe,
      mass = me,
      
      -- velocity space grid
      lower = {-4*vte, 0},
      upper = {4*vte, 12*me*vte^2/(4*B0)},
      cells = {10, 5},
      decompCuts = {1, 1},
      nDistFuncFrame = 16,   -- number of frames to output for distF
      
      -- initial conditions
      init = Plasma.MaxwellianProjection { 
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
	    return Te0
	 end,
      },

      -- collision parameters
      coll = Plasma.LBOCollisions {
	  collideWith = {'electron', 'ion'},
      },

      -- source parameters
      source = Plasma.MaxwellianProjection {
	 density = sourceDensity,
	 temperature = sourceTemperatureElc,
	 power = P_src/2,  -- source power evenly split bewtween ions and electrons
	 isSource = true,
      },

      -- boundary conditions
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},

      -- diagnostics
      diagnosticIntegratedMoments = { "intM0", "intM2"},
      diagnosticMoments = {"GkM0", "GkTemp", "GkUpar"}, 

      randomseed = randomseed,
   },

   -- gyrokinetic ions
   ion = Plasma.Species {
      evolve = true, -- evolve species?
      charge = qi,
      mass = mi,

      -- velocity space grid
      lower = {-6*vti, 0},
      upper = {6*vti, 12*mi*vti^2/(2*B0)},
      cells = {10, 5},
      decompCuts = {1, 1},
      nDistFuncFrame = 16,   -- number of frames to output for distF

      -- initial conditions
      init = Plasma.MaxwellianProjection { 
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
	    local x = xn[1]
	    return Ti0
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

      -- collision parameters
      coll = Plasma.LBOCollisions {
	 collideWith = {'ion', 'electron'},
      },

      -- source parameters
      source = Plasma.MaxwellianProjection {
	 density = sourceDensity,
	 temperature = sourceTemperatureIon,
	 power = P_src/2,  -- source power evenly split between ions and electrons
	 isSource = true,
      },

      -- boundary conditions
      bcx = {Plasma.Species.bcZeroFlux, Plasma.Species.bcZeroFlux},
      bcz = {Plasma.Species.bcSheath, Plasma.Species.bcSheath},
      
      -- diagnostics
      diagnosticMoments = {"GkM0", "GkTemp", "GkUpar"}, 
      diagnosticIntegratedMoments = { "intM0", "intM2"},

      randomseed = randomseed,
   },

   -- field solver
   field = Plasma.Field {
      -- dirichlet in x
      phiBcLeft = { T ="D", V = 0.0},
      phiBcRight = { T ="D", V = 0.0},
      -- periodic in y --
      -- no bc in z
      phiBcBack = { T ="N", V = 0.0},
      phiBcFront = { T ="N", V = 0.0},
      evolve = true, -- evolve fields?
   },

   -- magnetic geometry 
   funcField = Plasma.Geometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R0/x
      end,

      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
