-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
-- NSTX SOL test in 3x2v with 3x3v neutral species.
-- Run using Gkyl changeset [8a2b6a47522b].
----------------------------------------------------------------------

local Plasma = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi = require "Comm.Mpi"

-- scan parameters
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
nuFrac = 0.01
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
Lx = 50*rho_s
Ly = 100*rho_s
Lz = 8 -- [m]

-- Source parameters.
psolFac      = 0.25
P_SOL        = 5.4e6*nfac*Tfac*psolFac -- [W] 
P_src        = P_SOL*Ly*Lz/(2*math.pi*R*Lpol)
xSource      = R - 0.05 -- [m], source start coordinate
lambdaSource = 0.005 -- [m], characteristic length scale of density and temperature

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
   local sTfac = 4.0/psolFac
   if x < xSource + 3*lambdaSource then
      return sTfac*80*eV*Tfac
   else
      return sTfac*30*eV*Tfac
   end
end
sourceTemperatureIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local sTfac = 4.0/psolFac
   if x < xSource + 3*lambdaSource then
      return sTfac*80*eV*Tfac
   else
      return sTfac*30*eV*Tfac
   end
end
sourceDensityNeut = function (t,xn)
   return 8e-21*n0^2
end


-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = .5e-3,                        -- End time.
   nFrame      = 500,                         -- Number of output frames.
   lower       = {R - Lx/2-.02, -Ly/2, -Lz/2}, -- Configuration space lower left.
   upper       = {R + Lx/2, Ly/2, Lz/2},       -- Configuration space upper right.
   cells       = {30, 64, 32},                 -- Configuration space cells.
   basis       = "serendipity",                -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                            -- Polynomial order.
   timeStepper = "rk3",                        -- One of "rk2" or "rk3".
   cflFrac     = 0.5,
   restartFrameEvery = 0.001,

   -- Decomposition for configuration space.
   decompCuts = {6, 16, 8}, -- Cuts in each configuration direction.
   useShared = false,      -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {2},     -- Periodic in y only.

   -- Gyrokinetic electrons
   electron = Plasma.Species {
      evolve = true, -- Evolve species?
      charge = qe,
      mass  = me,

      -- Velocity space grid
      lower = {-4*vte, 0},
      upper = {4*vte, 12*me*vte^2/(2*B0)},
      cells = {12, 6},
      decompCuts = {1, 1},

      -- Initial conditions
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
                 local Ls              = Lz/4
                 local floor = 0.1
                 local effectiveSource = math.max(sourceDensity(t,{x,y,0}), floor)
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi*psolFac)
                 local nPeak           = 4.905066457001914e+23*4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
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
                    return 50*eV*Tfac
                 else 
                    return 20*eV*Tfac
                 end
              end,
              scaleWithSourcePower = false,
      },
      
      -- Source parameters
      source = Plasma.Source {
	 kind = "Maxwellian",
	 density = sourceDensity,
	 temperature = sourceTemperature,
	 power = P_src/2,
      },

      -- Collision parameters
      coll   = Plasma.LBOCollisions { 
         collideWith = {'electron', 'ion'},
         frequencies = {nuElc, nuElcIon},
         --nuFrac = nuFrac,
      },

      -- -- Neutral interactions
      ionization = Plasma.Ionization {
      	 collideWith  = {"neutral"},
      	 electrons    = "electron",
      	 neutrals     = "neutral",
      	 elemCharge   = eV, 
      	 elcMass      = me,
      	 plasma       = "H",
	 diagnostics = {"M0","intM0", "reactRate"},
      },

      -- Boundary conditions
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0"}},
             Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0"}}},

      -- Diagnostics
      diagnostics = { "M0", "M1", "M2", "Upar", "VtSq", "intM0", "intM1", "intM2"},
      nDistFuncFrame = 10,

      randomseed = randomseed,
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      evolve = true, -- Evolve species?
      charge = qi,

      mass   = mi,
      -- Velocity space grid.
      lower = {-4*vti, 0},
      upper = {4*vti, 12*mi*vti^2/(2*B0)},
      cells = {12, 6},
      decompCuts = {1, 1},

      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                 local x, y, z         = xn[1], xn[2], xn[3]
                 local Ls              = Lz/4
                 local floor = 0.1
                 local effectiveSource = math.max(sourceDensity(t,{x,y,0}), floor)
                 local c_ss            = math.sqrt(5/3*sourceTemperature(t,{x,y,0})/mi*psolFac)
                 local nPeak           = 4.905066457001914e+23*4*math.sqrt(5)/3/c_ss*Ls*effectiveSource/2
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
                    return 50*eV*Tfac
                 else 
                    return 20*eV*Tfac
                 end
              end,
              driftSpeed = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local Te
                 if x < xSource + 3*lambdaSource then 
                    Te = 50*eV
                 else 
                    Te = 20*eV
                 end
                 if math.abs(z) <= Lz/4 then
		    return z/(Lz/4)*math.sqrt(Te/mi)
                 else
		    return z/math.abs(z)*math.sqrt(Te/mi)
                 end
              end,
              scaleWithSourcePower = false,
      },

      -- Source parameters
      source = Plasma.Source {
	 kind = "Maxwellian",
	 density = sourceDensity,
	 temperature = sourceTemperatureIon,
	 power = P_src/2,
      },

      -- Collision parameters
      coll   = Plasma.LBOCollisions { 
         collideWith = {'ion', 'electron'},
         frequencies = {nuIon, nuIonElc},
         --nuFrac = nuFrac,
      },

      -- Neutral interactions
      ionization = Plasma.Ionization {
      	 collideWith  = {"neutral"},
      	 electrons    = "electron",
      	 neutrals     = "neutral",
      	 elemCharge   = eV,
      	 elcMass      = me,
      	 plasma       = "H",
      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"neutral"},
      	 ions = "ion",
      	 neutrals = "neutral",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "D",
      	 charge = qi,
      },

      -- Boundary conditions                                                                                                                                   
      bcx = {Plasma.ZeroFluxBC{}, Plasma.ZeroFluxBC{}},
      bcz = {Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0"}},
             Plasma.SheathBC{diagnostics={"M0","Upar","Energy","intM0"}}},

      -- Diagnostics                                                                                                                                           
      diagnostics = { "M0", "M1", "M2", "Upar", "VtSq", "intM0", "intM1", "intM2"},
      nDistFuncFrame = 10,

      randomseed = randomseed,
   },

   neutral = Plasma.Vlasov {
      evolve = true,
      charge = 0.0, 
      mass = mi,

      -- Velocity space grid
      lower = {-4.0*vti, -4.0*vti, -4.0*vti},
      upper = {4.0*vti, 4.0*vti, 4.0*vti},
      cells = {6, 6, 6},
      decompCuts = {1},

      -- Nonconservative positivity
      nonconPositivity = true,

      -- Initial conditions
      init = Plasma.VmMaxwellianProjection {
         density = function (t, xn)
   	    local x, y, z = xn[1], xn[2], xn[3]
            local n_n = n0
   	    local x0 = xSource + 3*lambdaSource
   	    local z0 = 0.1
   	    local flr = 0.01
   	    local fac = 1.0
   	    if z <= 0 then
   	       return fac*n_n*(sech((-Lz/2-z)/z0)^2 + flr)
   	    else	       
   	       return fac*n_n*(sech((Lz/2-z)/z0)^2 + flr)
   	    end
          end,
         driftSpeed = function (t, xn)
            return {0,0,0} --uPari
         end,
         temperature = function (t, xn)
            return 10*eV
         end,
      },

      -- Neutral interactions
      ionization = Plasma.Ionization {
      	 collideWith  = {"electron"},
      	 electrons    = "electron",
      	 neutrals     = "neutral",
      	 elemCharge   = eV, 
      	 elcMass      = me,
      	 plasma       = "H",
	 diagnostics = {"intM0"},
      },
      chargeExchange = Plasma.ChargeExchange {
      	 collideWith = {"ion"},
      	 ions = "ion",
      	 neutrals = "neutral",
      	 ionMass = mi,
      	 neutMass = mi,
      	 plasma = "D",
      	 charge = 0,
      },

      -- Source parameters
      source = Plasma.VmSource {
	 kind = "Maxwellian",
	 density = sourceDensityNeut,
	 driftSpeed = function (t, xn)
	    return {0,0,0}
	 end,
	 temperature = function (t, xn)
	    return 10*eV
	 end,
      },   

      -- Boundary conditions
      bcx = {Plasma.VmAbsorbBC{diagnostics={"M2"}}, Plasma.VmReflectBC{}},
      bcz = {Plasma.NeutralRecyclingBC{
                recycleTemp = 10*eV,  recycleIon  = "ion",
                recycleFrac = 0.95,
                diagnostics = {"M0","intM0","recycleCoef"},
             },
             Plasma.NeutralRecyclingBC{
                recycleTemp = 10*eV,  recycleIon  = "ion",
                recycleFrac = 0.95,  
                diagnostics = {"M0","intM0","recycleCoef"},
             }},

      -- Diagnostics=
      diagnostics = { "M0", "Udrift", "VtSq", "intM0", "intM1i", "intM2Flow", "intM2Thermal"},
      
      nDistFuncFrame = 10,
      
   },
   
   -- Field solver.
   field = Plasma.Field {
      -- Dirichlet in x.
      bcLowerPhi  = {{T = "D", V = 0.0}, {T = "P"}, {T ="N", V = 0.0}},
      bcUpperPhi  = {{T = "D", V = 0.0}, {T = "P"}, {T ="N", V = 0.0}},
      bcLowerApar = {{T ="D", V = 0.0}, {T = "P"}},
      bcUpperApar = {{T ="D", V = 0.0}, {T = "P"}},
      evolve     = true, -- Evolve fields?
      isElectromagnetic = false,
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
