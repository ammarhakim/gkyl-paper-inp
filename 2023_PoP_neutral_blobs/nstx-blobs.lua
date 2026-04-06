-- Gkyl ------------------------------------------------------------------------
--
-- NSTX-like gyrokinetic simulation with seeded blobs.
--
--
--------------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"
local Logger    = require "Lib.Logger"
local xsys      = require "xsys"

local log = Logger {
   logToFile = xsys.pickBool(logToFile, true)
}

-- Scan parameters.
nfac = 1
Tfac = 1

function sech(x)
   return 2*math.exp(x)/(math.exp(2*x)+1)
end
-- Universal constant parameters.
eps0         = Constants.EPSILON0
eV           = Constants.ELEMENTARY_CHARGE
qe           = -eV
qi           = eV
me           = Constants.ELECTRON_MASS

-- Plasma parameters.
mi           = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0          = 40*eV*Tfac
Ti0          = 40*eV*Tfac 
n0           = 7e18*nfac  -- [1/m^3]

-- Geometry and magnetic field.
B_axis       = 0.5   -- [T]
R0           = 0.85  -- [m]
a0           = 0.5   -- [m]
R            = R0 + a0
B0           = B_axis*(R0/R) -- [T]
Lpol         = 2.4 -- [m]

-- Parameters for collisions.
nuFrac = 0.1
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

log(string.format(" 1000/omega_ci = %e\n",2500/omega_ci))

-- Box size.
Lx = 100*rho_s
Ly = 100*rho_s
Lz = 8 -- [m]

-- Source parameters. (Not used for blob simulation, except xSource.)
P_SOL        = 5.4e6*nfac*Tfac -- [W] 
P_src        = P_SOL*Ly*Lz/(2*math.pi*R*Lpol)
xSource      = R - 0.05 -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature

-- Blob parameters
Lc = Lz
-- Blob initial conditions
a0b = (4*Lc^2/(rho_s*R))^(1/5)*rho_s
chi0 = xSource
aperp = 0.5*a0b
apar = Lz 
nb = 2*n0
Tb = 2*Te0
p_b = nb*2*Tb
delta_p = nb*2*Tb - n0*2*Te0

blobFunc = function (t,xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return math.exp(-math.log(2)*(((x - chi0)^2 + y^2)/(aperp^2) + z^2/apar^2))
end 

-- Multiply ion-gc density ICs by this factor to initialize a potential dipole
nIonFac = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return (-8.0*(math.log(2)^2)*y*rho_s^3/aperp^4*(2 -math.log(2)*((x-chi0)^2+y^2)/aperp^2) + 1)
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 10e-6,                         -- End time.
   nFrame      = 100,                            -- Number of output frames.
   lower       = {R - Lx/2, -Ly/2, -Lz/2}, -- Configuration space lower left.
   upper       = {R + Lx/2, Ly/2, Lz/2},       -- Configuration space upper right.
   cells       = {36, 36, 16},                    -- Configuration space cells.
   basis       = "serendipity",                -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                            -- Polynomial order.
   timeStepper = "rk3",                        -- One of "rk2" or "rk3".
   cflFrac     = 0.5,
   --restartFrameEvery = 0.01,
   calcIntQuantEvery = 1./10000.,

   groupDiagnostics = true,
   
   -- Decomposition for configuration space.
   decompCuts = {1, 1, 4}, -- Cuts in each configuration direction.
   useShared  = false,      -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {2},     -- Periodic in y only.

   -- Gyrokinetic electrons
   electron = Plasma.Species {
      evolve = true, -- Evolve species?
      charge = qe,
      mass   = me,
      -- Velocity space grid
      lower = {-4*vte, 0},
      upper = { 4*vte, 12*me*vte^2/(2*B0)},
      cells = {12, 6},
      decompCuts = {1, 1},

      -- Initial conditions
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
	    return n0 + nb*blobFunc(t, {x,y,z})
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3] 
	    return Te0 + Tb*blobFunc(t, {x,y,z})
	 end,
	 driftSpeed = function (t, xn)
	    local x, y, z = xn[1], xn[2], xn[3]
	    return c_s*(2*z/Lz)
	 end,
         scaleWithSourcePower = false,
      },
      
      -- Collision parameters
      coll = Plasma.LBOCollisions { 
         collideWith = {'electron', 'ion'},
         frequencies = {nuElc, nuElcIon},
      },

      -- -- Neutral interactions
      -- ionization = Plasma.Ionization {
      -- 	 collideWith = {"neutral"},  elemCharge = eV, 
      -- 	 electrons   = "electron",   elcMass    = me,
      -- 	 neutrals    = "neutral",    plasma     = "H",         
      -- },

      -- Boundary conditions
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{}, Plasma.SheathBC{}},

      -- Diagnostics
      diagnostics = {"M0", "Upar", "Temp"}, 
      nDistFuncFrame = 1,

      randomseed = randomseed,
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      evolve = true, -- Evolve species?
      charge = qi,
      mass   = mi,
      -- Velocity space grid.
      lower = {-4*vti, 0},
      upper = { 4*vti, 12*mi*vti^2/(2*B0)},
      cells = {12, 6},
      decompCuts = {1, 1},

      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
	    return n0 + nb*blobFunc(t, {x,y,z})*nIonFac(t, {x,y,z})
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3] 
	    return Ti0 + Tb*blobFunc(t, {x,y,z})
	 end,
	 driftSpeed = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
	    return c_s*(2*z/Lz)
         end,

         scaleWithSourcePower = false,
         -- driftSpeed = function (t, xn)
         --    local x, y, z = xn[1], xn[2], xn[3]
	 --    return c_s*(2*y/Lx - 1)
         -- end,
      },

      -- Source parameters.
      -- source = Plasma.Source {
      --    kind        = "Maxwellian",
      --    density     = sourceDensity,
      --    temperature = sourceTemperature,
      --    power       = P_src/2,
      -- },

      -- Collision parameters.
      coll = Plasma.LBOCollisions { 
         collideWith = {'ion', 'electron'},
         frequencies = {nuIon, nuIonElc},
      },

      -- -- Neutral interactions.
      -- ionization = Plasma.Ionization {
      -- 	 collideWith  = {"neutral"},  elemCharge   = eV,
      -- 	 electrons    = "electron",   elcMass      = me,
      -- 	 neutrals     = "neutral",    plasma       = "H",
      -- },
      -- chargeExchange = Plasma.ChargeExchange {
      -- 	 collideWith = {"neutral"},  neutMass = mi,
      -- 	 ions        = "ion",        plasma   = "H",
      -- 	 neutrals    = "neutral",    charge   = qi,
      -- 	 ionMass     = mi,
      -- },

      -- Boundary conditions.
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{}, Plasma.SheathBC{}},

      -- Diagnostics.
      diagnostics = {"M0", "Upar", "Temp"}, 
      nDistFuncFrame = 1,

      randomseed = randomseed,
   },

--    neutral = Plasma.Vlasov {
--       evolve = true,
--       charge = 0.0, 
--       mass   = mi,
--       -- Velocity space grid.
--       lower = {-4.0*vti, -4.0*vti, -4.0*vti},
--       upper = { 4.0*vti,  4.0*vti,  4.0*vti},
--       cells = {4, 4, 4},
--       decompCuts = {1},

--       -- Initial conditions.
--       init = Plasma.VmMaxwellianProjection {
--          density = function (t, xn)
--    	    local x, y, z = xn[1], xn[2], xn[3]
--             local n_n = 0.1*n0
--    	    local x0  = xSource + 3*lambdaSource
--    	    local z0  = 0.1
--    	    local flr = 0.01
--    	    local fac = 1.0
--    	    if x < x0 then
--    	       fac = math.exp((x-x0)/0.01)
--    	    end
--    	    if z <= 0 then
--    	       return fac*n_n*(sech((-Lz/2-z)/z0)^2 + flr)
--    	    else	       
--    	       return fac*n_n*(sech((Lz/2-z)/z0)^2 + flr)
--    	    end
--           end,
--          driftSpeed = function (t, xn)
--             return {0,0,0} --uPari
--          end,
--          temperature = function (t, xn)
--             return 2*eV
--          end,
--       },

--       -- Neutral interactions.
--       ionization = Plasma.Ionization {
--       	 collideWith = {"electron"},  elemCharge = eV, 
--       	 electrons   = "electron",    elcMass    = me,
--       	 neutrals    = "neutral",     plasma     = "H",         
--       },
--       chargeExchange = Plasma.ChargeExchange {
--       	 collideWith = {"ion"},    neutMass = mi,
--       	 ions        = "ion",      plasma   = "H",
--       	 neutrals    = "neutral",  charge   = 0,
--       	 ionMass     = mi,
--       },

--       -- Source parameters.
--       -- source = Plasma.VmMaxwellianProjection{
--       source = Plasma.VmSource {
--          kind        = "Maxwellian",
--          density     = sourceDensityNeut,
--          temperature = 2.*eV,
--       },

--       -- Boundary conditions.
--       bcx = {Plasma.AbsorbBC{}, Plasma.AbsorbBC{}},
--       bcz = {Plasma.ReflectBC{}, Plasma.ReflectBC{}},

--       -- Diagnostics.
--       diagnostics = { "M0", "Udrift", "VtSq"},
--       -- diagnostics = { "M0", "Udrift", "VtSq", "intM0", "intM1i", "intM2Flow", "intM2Thermal"},
--    },
   
   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
      isElectromagnetic = false,
      -- Dirichlet in x, periodic in y. Potential phi has homogeneous Neumann
      -- BC for the smoothing operation that enforces continuity in z.
      bcLowerPhi = {{T = "D", V = 0.0}, {T = "P"}, {T ="N", V = 0.0}},
      bcUpperPhi = {{T = "D", V = 0.0}, {T = "P"}, {T ="N", V = 0.0}},
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
-- Run application.
plasmaApp:run()
