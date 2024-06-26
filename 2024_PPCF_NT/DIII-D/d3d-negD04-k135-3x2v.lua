-- Gkyl ------------------------------------------------------------------------
--
-- An 3x2v open field line simulation with Miller geometry using
-- parameters similar to those near the LCFS of LTX's shot 103795
-- at t=469.11 ms.
--
-- The plasma parameters are taken from Elizabeth Perez's old Li simulation.
-- ~/gkylsoft/openmpi/bin/mpirun -n 4 ~/gkylsoft/gkyl/bin/gkyl *lua
--
--------------------------------------------------------------------------------
local Plasma    = (require "App.PlasmaOnCartGrid").Gyrokinetic()
local Constants = require "Lib.Constants"
local Mpi       = require "Comm.Mpi"
local xsys      = require "xsys"
local Logger    = require "Lib.Logger"
local math      = require("sci.math").generic
local root      = require("sci.root")
local quad      = require("sci.quad")
local diff      = require("sci.diff-recursive")
local df        = diff.df

local log = Logger {
   logToFile = xsys.pickBool(logToFile, true)
}

-- Use the following to disable embedding of the input file in our data files.
-- Otherwise that tries to store a very large string in an input file which
-- ADIOS or some filesystems don't robustly do and may cause I/O errors.
GKYL_EMBED_INP = false

-- Universal constant parameters.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
qe, qi   = -eV, eV
me, mp   = Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters. Chosen based on the value of a cubic sline
-- between the last TS data inside the LCFS and the probe data in
-- in the far SOL, near R=0.475 m.
AMU = 2.01410177811
mi  = mp*AMU  -- Deuterium ions.
Te0 = 100*eV 
Ti0 = 100*eV 
n0  = 2.0e19     -- [1/m^3]

-- Geometry and magnetic field.
--R_axisTrue = 0.388252435        -- [m]
Rdim       = 1.7  -- [m]
Zdim       = 3.2  -- [m]
Z_axis     = 0.013055028 -- [m]
R_axisTrue = 1.6486461 --0.95*R_axisTrue    -- Change R_axis to fit geometry better.
R_axis     = 1.6
B_axis     = 2*R_axisTrue/R_axis   -- [T]
R_LCFSmid  = 2.17 --2.17885               -- Major radius of the LCFS at the outboard midplane [m].
Rmid_min   = R_LCFSmid - 0.1          -- Minimum midplane major radius of simulation box [m].
Rmid_max   = R_LCFSmid + 0.05 --2.32               -- Maximum midplane major radius of simulation box [m].
R0         = 0.5*(Rmid_min+Rmid_max)  -- Major radius of the simulation box [m].
a_mid      = R_LCFSmid-R_axis   -- Minor radius at outboard midplane [m].
r0         = R0-R_axis          -- Minor radius of the simulation box [m].
B0         = B_axis*(R_axis/R0) -- Magnetic field magnitude in the simulation box [T].

-- What are these values??
qSep       = 5.22              -- Safety factor at the separatrix.
sSep       = 1.27976219        -- Magnetic shear at the separatrix.
kappa      = 1.35 -- 1.33              -- Elongation (=1 for no elongation).
delta      = -0.4 --0.4             -- Triangularity (=0 for no triangularity).

function r_x(x) return x+a_mid-0.1 end ---0.1 end   -- Minor radius given x.
--function qprofile(r) return qSep/(1.-sSep*((r-a_mid)/a_mid)) end   -- Magnetic safety profile.
function qprofile(r)
   local a = {49.46395467479657, -260.79513158768754, 458.42618139184754, -267.63441353752336}
   return a[1]*(r+R_axis+0.15)^3 + a[2]*(r+R_axis+0.15)^2 + a[3]*(r+R_axis+0.15) + a[4]
end

nuFrac = 0.1
-- Electron collision freq.
logLambdaElc = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Te0/eV)
nuElc = nuFrac*logLambdaElc*eV^4*n0/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(me)*(Te0)^(3/2))
-- Ion collision freq.
logLambdaIon = 6.6 - 0.5*math.log(n0/1e20) + 1.5*math.log(Ti0/eV)
nuIon = nuFrac*logLambdaIon*eV^4*n0/(12*math.pi^(3/2)*eps0^2*math.sqrt(mi)*(Ti0)^(3/2))
-- Electron-ion and ion-electron collision frequencies.
nuElcIon = nuElc*math.sqrt(2)
nuIonElc = nuElcIon/(mi/me)

-- Derived parameters
vti, vte = math.sqrt(Ti0/mi), math.sqrt(Te0/me)
c_s      = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
rho_s    = c_s/omega_ci

-- Box size.
Lx = Rmid_max-Rmid_min
xMin, xMax = 0., Lx
rMin, rMax = Rmid_min-R_axis, Rmid_max-R_axis
q0 = qprofile(r_x(0.5*(xMin+xMax)))    -- Magnetic safety factor in the center of domain.
--ntoroidal = 6
Ly = 150*rho_s --2*math.pi*(r0/q0)/ntoroidal --50*rho_s
Lz = 2.*math.pi+1e-8
x_LCFS = R_LCFSmid - Rmid_min
ntoroidal = 2*math.pi*r0/q0/Ly
log(string.format(" Lx = %f\n",Lx))
log(string.format(" Ly = %f\n",Ly))
log(string.format(" Lz = %f\n",Lz))
log(string.format(" x_LCFS = %f\n",x_LCFS))
log(string.format(" Lx/rho_s = %f\n",Lx/rho_s))

epsilon0  = r0/R0              -- Inverse aspect ratio in the center of the domain.
nuStarElc = nuElc*q0*R0/(vte*(epsilon0^(3./2.)))
nuStarIon = nuIon*q0*R0/(vti*(epsilon0^(3./2.)))
log(string.format(" nuStarElc = %g\n", nuStarElc))
log(string.format(" nuStarIon = %g\n", nuStarIon))

-- Functions needed for (x,y,z) -> (X,Y,Z) mapping in Miller geometry.
local function R(r, theta) return R_axis + r*math.cos(theta + math.asin(delta)*math.sin(theta)) end
local function Z(r, theta) return Z_axis + kappa*r*math.sin(theta) end
local function Bphi(R) return B0*R0/R end
local function Jr(r, theta)
   return R(r,theta)*(df(R,1)(r,theta)*df(Z,2)(r,theta) - df(Z,1)(r,theta)*df(R,2)(r,theta))
end
local function dPsidr(r, theta)
   local function integrand(t) return Jr(r,t)/R(r,t)^2 end
   local integral
   integral, _ = quad.dblexp(integrand, 0, 2*math.pi, 1e-10)
   return B0*R_axis/(2*math.pi*qprofile(r))*integral
end
local function J(r, theta) return Jr(r, theta)/dPsidr(r, theta) end
local function alpha(r, theta, phi)
   local function integrand(t) return Jr(r,t)/R(r,t)^2 end
   local integral
   local t = theta
   while diff.lt(t, -math.pi) do t = t+2*math.pi end
   while diff.lt( math.pi, t) do t = t-2*math.pi end
   if diff.lt(0, t) then
      integral, _ =  quad.dblexp(integrand, 0, t, 1e-10)/dPsidr(r,theta)
   else
      integral, _ = -quad.dblexp(integrand, t, 0, 1e-10)/dPsidr(r,theta)
   end
   return phi - B0*R_axis*integral
end
local function gradr(r, theta)
   return R(r,theta)/Jr(r,theta)*math.sqrt(df(R,2)(r,theta)^2 + df(Z,2)(r,theta)^2)
end

-- Source parameters.
-- TRANSP estimated P_OH=1.25e5 W for this shot.
P_SOL     = 1.5e6           -- Power into the whole SOL, from experimental heating power [W].
P_src     = P_SOL/ntoroidal -- Amount of power into flux tube [W].
n0_src    = 6.6e23          -- Amplitude of density source [m^{-3} s^{-1}].
x_src     = xMin --+0.3*Lx     -- Source start coordinate [m].
sigma_src = 0.03*Lx         -- Characteristic length scale of source [m].
T_srce    = 2*Te0           -- Electron source temperature [J].
T_srci    = 2*Ti0           -- Ion source temperature [J].
x_srcGB   = xMin
sig_srcGB = 10*rho_s
sigma_init = 0.1*Lx
SN0 = 9e22 --7.786691893452529e+22, 1.0646985811428132e+23
SN_GB = 8.092675420182799e+21*1.1
bFac = 1.2

-- Source profiles.
initDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local srcFloor = 1e-2
   return 1e19*(math.exp(-((x-x_src)^2)/(2.*(sigma_init^2)))+srcFloor)
end
sourceDensity = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local srcFloor = 1e-2
   if x < x_src then srcFloor = 1e-2 end -- Higher floor to left of source peak.
   --if math.abs(z) < Lz/4 then
   return SN0*(math.exp(-((x-x_src)^2)/(2.*(sigma_src^2)))+srcFloor)
   -- else
   --    return 1e-40
   -- end
end
sourceDensElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local srcFloor = 1e-2
   -- local SN = 1.639951372314892e+23 -- 8.506198416096726e+21
   -- Source with poloidal dependence
   return math.exp(-(x-x_srcGB)^2/(2.*sig_srcGB^2))
      *math.max(SN_GB*math.sin(z)*math.exp(-math.abs(z)^1.5/(2*bFac^2)),0.)
end
sourceDensIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local srcFloor = 1e-2
   -- local SN = -1.639951372314892e+23 -- -8.506198416096726e+21
   -- Source with poloidal dependence
   return math.exp(-(x-x_srcGB)^2/(2.*sig_srcGB^2))
      *math.max(-SN_GB*math.sin(z)*math.exp(-math.abs(z)^1.5/(2*bFac^2)),0.)
end
sourceTempIon = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < x_src + 3*sigma_src then
      return T_srci
   else
      return T_srci*3./8.
   end
end
sourceTempElc = function (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   if x < x_src + 3*sigma_src then
      return T_srce
   else
      return T_srce*3./8.
   end
end
sourceTempIonGB = function (t, xn)
   return 100*eV
end
sourceTempElcGB = function (t, xn)
   return 100*eV
end

-- Initialize a random seed for initial conditions
-- will be used for both ions and electrons.
randomseed = 100000*Mpi.Comm_rank(Mpi.COMM_WORLD)+63--os.time()

local bcShiftFunc = function(t,xn)
   local x, y, z = xn[1], xn[2], xn[3]
   local r = r_x(x)
   return r0/q0*qprofile(r)*Lz
end

local function stop(tol)
   return function(x, y, xl, xu, yl, yu)
      if diff.lt(math.abs(y), tol) then return true
      else return false
      end
   end
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = 1.e-3,                 -- End time.
   nFrame = 1000,                  -- Number of output frames.
   lower  = {xMin,-Ly/2,-Lz/2},                  -- Configuration space lower left.
   upper  = {xMax, Ly/2, Lz/2},                  -- Configuration space upper right.
   cells  = {96, 96, 16},                     -- Configuration space cells.
   mapc2p = function(xc)                   -- Transformation from computational to physical coordinates.
      local x, y, z = xc[1], xc[2], xc[3]
      local r = r_x(x)
      -- Map to cylindrical (R, Z, phi) coordinates.
      local R   = R(r, z)
      local Z   = Z(r, z) --kappa*r*math.sin(z)
      local phi = -q0/r0*y - alpha(r, z, 0)
      -- Map to Cartesian (X, Y, Z) coordinates.
      local X = R*math.cos(phi)
      local Y = R*math.sin(phi)
      return X, Y, Z
   end,
   basis       = "serendipity",            -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                        -- Polynomial order.
   timeStepper = "rk3",                    -- One of "rk2" or "rk3".
   cflFrac     = 0.6,
   restartFrameEvery = .01,
   calcIntQuantEvery = 1./(1000*10.),   -- Aim at 10x more frequently than frames.
   groupDiagnostics  = true,
   
   periodicDirs = {2},     -- Periodic in y only.

   decompCuts = {12,12,4},    -- MPI subdomains/processes.

   -- Gyrokinetic electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      lower = {-4.*vte, 0},
      upper = { 4.*vte, me*((4.*vte)^2)/(2.*B0)},
      cells = {12, 6},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z         = xn[1], xn[2], xn[3]
            -- local Ls              = Lz/4.
            -- local floor           = 0.1
            -- local effectiveSource = math.max(initDensity(t,{x,y,0}), floor)
            -- --local c_ss            = math.sqrt((5./3.)*sourceTempElc(t,{x,y,0})/mi)
            -- --local nPeak           = (4.*math.sqrt(5.)/3.)*(Ls/c_ss)*effectiveSource/2.
            -- local perturb         = 0 
            -- if math.abs(z) <= Ls then
            --    return effectiveSource*(1.+math.sqrt(1.-(z/Ls)^2))/2.*(1.+perturb)
            -- else
            --    return (effectiveSource/2.)*(1.+perturb)
            -- end
	    return 1e19*(0.5*(1.+math.tanh(2.*(2.-25.*x)))+0.01)
         end,
         temperature = function (t, xn)
            local x = xn[1]
            -- if (x < x_src + 3*sigma_src) then
            --    return Te0*1.25
            -- else
            --    return Te0/2.
            -- end
	    return Te0*((1/3)*(2.+math.tanh(2.*(2.-25.*x)))+0.01)
         end,
         --scaleWithSourcePower = true,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc','ion'},
         frequencies = {nuElc,nuElcIon},
      },
      source = Plasma.Source {
         kind        = "Maxwellian",  density     = sourceDensity,
         --power = P_src/2,
	 temperature = sourceTempElc,
         diagnostics = {"M0","intM0"},
      },
      sourceGB = Plasma.Source {
         kind        = "Maxwellian",  density     = sourceDensElc,
      	 temperature = sourceTempElcGB,
         diagnostics = {"M0","intM0"},
      },
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "intM0", "intM1", "intEnergy",}, 
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      bcz = {Plasma.TokamakEdgeBC{xLCFS=x_LCFS, shiftFunction=bcShiftFunc,
                                  diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.TokamakEdgeBC{xLCFS=x_LCFS, shiftFunction=bcShiftFunc,
                                  diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      -- bcz = {Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
      --        Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
   },

   -- Gyrokinetic ions
   ion = Plasma.Species {
      charge = qi,  mass = mi,
      -- Velocity space grid.
      lower = {-4.*vti, 0},
      upper = { 4.*vti, mi*(4.*vti^2)/(2.*B0)},
      cells = {12, 6},
      -- Initial conditions.
      init = Plasma.MaxwellianProjection {
         density = function (t, xn)
            local x, y, z         = xn[1], xn[2], xn[3]
            -- local Ls              = Lz/4.
            -- local floor           = 0.1
            -- local effectiveSource = math.max(initDensity(t,{x,y,0}), floor)
            -- --local c_ss            = math.sqrt((5./3.)*sourceTempElc(t,{x,y,0})/mi)
            -- --local nPeak           = (4.*math.sqrt(5.)/3.)*(Ls/c_ss)*effectiveSource/2.
            -- local perturb         = 0
            -- if math.abs(z) <= Ls then
            --    return effectiveSource*(1.+math.sqrt(1.-(z/Ls)^2))/2.*(1.+perturb)
            -- else
            --    return (effectiveSource/2.)*(1.+perturb)
            -- end
	    return 1e19*(0.5*(1.+math.tanh(2.*(2.-25.*x)))+0.01)
         end,
         temperature = function (t, xn)
            local x = xn[1]
            -- if (x < x_src + 3*sigma_src) then
            --    return Ti0*1.25
            -- else
            --    return Ti0/2.
            -- end
	    return Ti0*((1/3)*(2.+math.tanh(2.*(2.-25.*x)))+0.01)
         end,
         --scaleWithSourcePower = true,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion','elc'},
         frequencies = {nuIon,nuIonElc},
      },
      source = Plasma.Source {
         kind        = "Maxwellian",  density     = sourceDensity,
         --power = P_src/2,
	 temperature = sourceTempIon,
         diagnostics = {"M0","intM0"},
      },
      sourceGB = Plasma.Source {
         kind        = "Maxwellian",  density     = sourceDensIon,
      	 temperature = sourceTempIonGB,
         diagnostics = {"M0","intM0"},
      },
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "intM0", "intM1", "intEnergy"}, 
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      bcz = {Plasma.TokamakEdgeBC{xLCFS=x_LCFS, shiftFunction=bcShiftFunc,
                                  diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.TokamakEdgeBC{xLCFS=x_LCFS, shiftFunction=bcShiftFunc,
                                  diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      -- bcz = {Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
      --        Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
   },

   -- Field solver.
   field = Plasma.Field {
      evolve = true, -- Evolve fields?
      -- Dirichlet in x, periodic in y. Potential phi has homogeneous Neumann
      -- BC for the smoothing operation that enforces continuity in z.
      bcLowerPhi = {{T = "D", V = 0.0}, {T = "P"}, {T = "N", V = 0.0}},
      bcUpperPhi = {{T = "D", V = 0.0}, {T = "P"}, {T = "N", V = 0.0}},
      isElectromagnetic = false,
   },

   -- Magnetic geometry.
   externalField = Plasma.Geometry {
      -- Background magnetic field.
      bmag = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         local r = r_x(x)
         local Bt = Bphi(R(r,z))
         local Bp = dPsidr(r,z)/R(r,z)*gradr(r,z)
         return math.sqrt(Bt^2 + Bp^2)
      end,

      -- Geometry is not time-dependent.
      evolve = false,
   },
}
-- Run application.
plasmaApp:run()
