-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"

-- physical parameters
gasGamma = 5./3.
epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

Te_Ti = 1.0 -- ratio of electron to ion temperature
n0 = 1.0 -- initial number density

elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge
ionMass = 100.0 -- ion mass
ionCharge = 1.0 -- ion charge

vAe = 0.2 --electron Alfven speed (how non-relativistic the system is)
B0 = vAe*math.sqrt(mu0*n0*elcMass)
beta = 1.0 --proton plasma beta = 2 mu0 ni Ti/B0^2
vtElc = vAe*math.sqrt(beta) --electron thermal velocity

elcTemp = vtElc^2/2.0
ionTemp = elcTemp/Te_Ti

-- ion velocities
vAi = vAe/math.sqrt(ionMass)
vtIon = vtElc/math.sqrt(ionMass) --Ti/Te = 1.0

-- cyclotron frequencies
omegaCe = ionCharge*B0/elcMass
omegaCi = ionCharge*B0/ionMass

-- inertial length
de = vAe/omegaCe
di = vAi/omegaCi

-- gyroradii
rhoi = vtIon/omegaCi
rhoe = vtElc/omegaCe

-- OT initial conditions
elong = 5.
u0x = vAi/elong
u0y = vAi/elong
B0x = B0/elong
B0y = B0/elong
ni0 = n0 -- Guess for initial value of ni0


-- domain size and simulation time
Lx = 100.0*math.pi*rhoi
Ly = Lx
Lz = elong*Lx
nx = 448
ny = 448
nz = 448

ncx = 56
ncy = 8
ncz = 8

tEnd = 1500.0/omegaCi
nOut = 75

momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nOut,
   lower = {0.0, 0.0, 0.0}, -- configuration space lower left
   upper = {Lx, Ly, Lz}, -- configuration space upper right
   cells = {nx, ny, nz}, -- configuration space cells
   timeStepper = "fvDimSplit",
   cfl = 1.0,
   -- decomposition for configuration space
   decompCuts = {ncx, ncy, ncz}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2,3}, -- periodic directions

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      -- initial conditions
      init = function (t, xn)
	 local x, y, z = xn[1], xn[2], xn[3]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge
	 local ne = n0

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

	 -- initial current computed from curl(B) = mu0 J
         -- local Jz = (B0y*(_4pi/Lx)*cos(_4pi*x/Lx) + B0x*(_2pi/Ly)*cos(_2pi*y/Ly)) / mu0

         -- local vdrift_x = -u0x*sin(_2pi*y/Ly)
         -- local vdrift_y = u0y*sin(_2pi*x/Lx)
         -- local vdrift_z = -Jz / qi
         
 	 local Jx = (_2pi/Lz)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) + cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jx = Jx + (_2pi/Lz)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) - cos(_4pi*x/Lx + _2pi*z/Lz)) / mu0
	 local Jy = (_2pi/Lz)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0
	 local Jz = (_2pi/Lx)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) - cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jz = Jz + (_4pi/Lx)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) + cos(_4pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jz = Jz + (_2pi/Ly)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0

         local vdrift_x = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz) - Jx / (qi*ne)
         local vdrift_y = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz)) - Jy / (qi*ne)
         vdrift_y = vdrift_y + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz))
         local vdrift_z = -Jz / (qi*ne)

	 local rhoe = ne*elcMass
         local exmom = vdrift_x*rhoe
	 local eymom = vdrift_y*rhoe
         local ezmom = vdrift_z*rhoe
	 local ere = ne*elcTemp/(gasGamma-1) + 0.5*exmom*exmom/rhoe + 0.5*eymom*eymom/rhoe + 0.5*ezmom*ezmom/rhoe
	 
	 return rhoe, exmom, eymom, ezmom, ere
      end,
      evolve = true, -- evolve species?
   },

   -- ions
   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,

      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      -- initial conditions
      init = function (t, xn)
	 local x, y, z = xn[1], xn[2], xn[3]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         -- no correction to ion density to initialize electric field with minimal transient
         local ni = n0

         -- local vdrift_x = -u0x*sin(_2pi*y/Ly)
         -- local vdrift_y = u0y*sin(_2pi*x/Lx)
         -- local vdrift_z = 0.0

       	 local vdrift_x = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz)
         local vdrift_y = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz))
         vdrift_y = vdrift_y + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz))
	 local vdrift_z = 0.0

	 local rhoi = ni*ionMass
         local ixmom = vdrift_x*rhoi
	 local iymom = vdrift_y*rhoi
         local izmom = vdrift_z*rhoi
	 local eri = ni*ionTemp/(gasGamma-1) + 0.5*ixmom*ixmom/rhoi + 0.5*iymom*iymom/rhoi + 0.5*izmom*izmom/rhoi
	 
	 return rhoi, ixmom, iymom, izmom, eri

      end,
      evolve = true, -- evolve species?
   },

   field = Moments.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
	 local x, y, z = xn[1], xn[2], xn[3]

         local cos = math.cos
         local sin = math.sin

         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

         local ne = n0
         local ni = n0

 	 local Jx = (_2pi/Lz)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) + cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jx = Jx + (_2pi/Lz)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) - cos(_4pi*x/Lx + _2pi*z/Lz)) / mu0
	 local Jy = (_2pi/Lz)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0
	 local Jz = (_2pi/Lx)*0.5*B0y*(cos(_2pi*x/Lx - _2pi*z/Lz) - cos(_2pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jz = Jz + (_4pi/Lx)*0.5*B0y*(cos(_4pi*x/Lx - _2pi*z/Lz) + cos(_4pi*x/Lx + _2pi*z/Lz) ) / mu0
	 Jz = Jz + (_2pi/Ly)*B0x*cos(_2pi*y/Ly - _2pi*z/Lz) / mu0


         local Bx = -B0x*sin(_2pi*y/Ly - _2pi*z/Lz)
	 local By = 0.5*B0y*(sin(_2pi*x/Lx - _2pi*z/Lz) - sin(_2pi*x/Lx + _2pi*z/Lz))
         By = By + 0.5*B0y*(sin(_4pi*x/Lx - _2pi*z/Lz) + sin(_4pi*x/Lx + _2pi*z/Lz))
         local Bz = B0

         -- Assumes qi = abs(qe)
	 local u_xe = -u0x*sin(_2pi*y/Ly - _2pi*z/Lz) - Jx / (qi*ne)
         local u_ye = 0.5*u0y*(sin(_2pi*x/Lx + _2pi*z/Lz) + sin(_2pi*x/Lx - _2pi*z/Lz)) - Jy / (qi*ne)
         u_ye = u_ye + 0.5*u0y*(sin(_4pi*x/Lx - _2pi*z/Lz) - sin(_4pi*x/Lx + _2pi*z/Lz))
         local u_ze = -Jz / (qi*ne)

         -- E = - v_e x B ~  (J - u) x B
         local Ex = - (u_ye*Bz - u_ze*By)
         local Ey = - (u_ze*Bx - u_xe*Bz)
         local Ez = - (u_xe*By - u_ye*Bx)

         return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- evolve field?
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
   },   

}
-- run application
momentApp:run()
