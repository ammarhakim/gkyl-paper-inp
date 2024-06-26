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

-- Universal constant parameters.
eps0, eV = Constants.EPSILON0, Constants.ELEMENTARY_CHARGE
qe, qi   = -eV, eV
me, mp   = Constants.ELECTRON_MASS, Constants.PROTON_MASS

-- Plasma parameters. 
AMU = 2.01410177811
mi  = mp*AMU  -- Deuterium ions.
Te0 = 40*eV 
Ti0 = 40*eV 
n0  = 7.0e18     -- [1/m^3]

-- Geometry and magnetic field.
--R_axisTrue = 0.388252435        -- [m]
Rdim       = 1.7  -- [m]
Zdim       = 3.2  -- [m]
Z_axis     = 0.013055028 -- [m]
--R_axisTrue = 1.6486461 --0.95*R_axisTrue    -- Change R_axis to fit geometry better.
R_axis     = 1.65
B_axis     = 2 --*R_axisTrue/R_axis   -- [T]
R_LCFSmid  = 2.17 --2.17885               -- Major radius of the LCFS at the outboard midplane [m].
Rmid_min   = R_LCFSmid          -- Minimum midplane major radius of simulation box [m].
Rmid_max   = 2.35 --2.32               -- Maximum midplane major radius of simulation box [m].
R0         = 0.5*(Rmid_min+Rmid_max)  -- Major radius of the simulation box [m].
a_mid      = R_LCFSmid-R_axis   -- Minor radius at outboard midplane [m].
r0         = R0-R_axis          -- Minor radius of the simulation box [m].
B0         = B_axis*(R_axis/R0) -- Magnetic field magnitude in the simulation box [T].

-- What are these values??
qAxis      = 0.7
qSep       = 3.0              -- Safety factor at the separatrix.
--sSep       = 1.27976219        -- Magnetic shear at the separatrix.
kappa      = 1.0 -- 1.33              -- Elongation (=1 for no elongation).
delta      = 0.0 --0.4             -- Triangularity (=0 for no triangularity).

function r_x(x) return x+a_mid end   -- Minor radius given x.
function qprofile(r) return qAxis + (qSep - qAxis)*(r/a_mid)^2 end   -- Magnetic safety profile.

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
Ly = 300*rho_s
Lz = 0.9*2.*math.pi+1e-8
log(string.format(" Lx = %f\n",Lx))
log(string.format(" Ly = %f\n",Ly))
log(string.format(" Lz = %f\n",Lz))

q0 = qprofile(r_x(0.5*(xMin+xMax)))    -- Magnetic safety factor in the center of domain.
ntoroidal = 2*math.pi*r0/q0/Ly
log(string.format("ntoroidal = %g\n", ntoroidal))

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

-- Blob parameters
Lc = 50
-- Blob initial conditions
a0b = (4*Lc^2/(rho_s*Rmid_min))^(1/5)*rho_s
chi0 = Lx/3
aperp = 1.0*a0b
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

local function stop(tol)
   return function(x, y, xl, xu, yl, yu)
      if diff.lt(math.abs(y), tol) then return true
      else return false
      end
   end
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd   = 10e-6, --1000/omega_ci,                 -- End time.
   nFrame = 100,                  -- Number of output frames.
   lower  = {xMin,-Ly/2,-Lz/2},                  -- Configuration space lower left.
   upper  = {xMax, Ly/2,Lz/2},                  -- Configuration space upper right.
   cells  = {72, 72, 16},                     -- Configuration space cells.
   mapc2p = function(xc)                   -- Transformation from computational to physical coordinates.
      local x, y, z = xc[1], xc[2], xc[3]
      local r = r_x(x)
      -- Map to cylindrical (R, Z, phi) coordinates.
      local R   = R(r, z)
      local Z   = Z(r, z) --kappa*r*math.sin(z)
      local phi = qprofile(r)/r0*y + alpha(r, z, 0)
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

   decompCuts = {12,16,4},    -- MPI subdomains/processes.

   -- Gyrokinetic electrons.
   elc = Plasma.Species {
      charge = qe,  mass = me,
      lower = {-4.*vte, 0},
      upper = { 4.*vte, me*((4.*vte)^2)/(2.*B0)},
      cells = {12, 6},
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
      coll = Plasma.LBOCollisions {
         collideWith = {'elc','ion'},
         frequencies = {nuElc,nuElcIon},
      },
      -- No sources
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "intM0", "intM1", "intEnergy",}, 
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      bcz = {Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
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
            local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
	            return n0 + nb*nIonFac(t, {x,y,z})*blobFunc(t, {x,y,z})
         end,
         temperature = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3] 
	            return Ti0 + Tb*blobFunc(t, {x,y,z})
		    	    end,
         scaleWithSourcePower = false,
         driftSpeed = function (t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
	    return c_s*(2*z/Lz)
         end,
      },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion','elc'},
         frequencies = {nuIon,nuIonElc},
      },
      -- No sources
      evolve = true, -- Evolve species?
      diagnostics = {"M0", "Upar", "Temp", "Tpar", "Tperp", "intM0", "intM1", "intEnergy"}, 
      nDistFuncFrame = 10,
      randomseed = randomseed,
      bcx = {Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.AbsorbBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
      bcz = {Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}},
             Plasma.SheathBC{diagnostics={"M0","Energy","intM0","intM1","intKE","intEnergy"}}},
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
