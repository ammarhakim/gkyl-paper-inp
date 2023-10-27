-- Gkyl ------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = true
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

epsilon0 = 1.0 -- permittivity of free space
mu0 = 1.0 -- pemiability of free space
lightSpeed = 1/math.sqrt(mu0*epsilon0) -- speed of light

Te_Ti = 0.2 -- ratio of electron to ion temperature
n0 = 1.0 -- initial number density in layer
n0_ninf = 5. -- upstream density, ninf
vAe = 1.0/4.0 -- Based on in-plane field only!
beta = 5./6.
-- plasma beta for electrons
beta_electron = beta*Te_Ti

elcMass = 1.0 -- electron mass
elcCharge = -1.0 -- electron charge
ionMass = 25.0 -- ion mass
ionCharge = 1.0 -- ion charge

B0 = vAe*math.sqrt(mu0*n0*elcMass)
Bg = 0.1*B0 
BTot = math.sqrt(B0^2 + Bg^2) 
vAeTot = BTot/math.sqrt(mu0*n0*elcMass) 

vtElc = vAeTot*math.sqrt(beta*Te_Ti)

elcTemp = vtElc^2/2.0
ionTemp = elcTemp/Te_Ti

-- ion velocities
vAi = vAe/math.sqrt(ionMass)
vtIon = vtElc/math.sqrt(ionMass*Te_Ti)

-- plasma and cyclotron frequencies
wpe = math.sqrt(ionCharge^2*n0/(epsilon0*elcMass))
wpi = math.sqrt(ionCharge^2*n0/(epsilon0*ionMass))
omegaCe = ionCharge*B0/elcMass
omegaCi = ionCharge*B0/ionMass
omegaCeTot = ionCharge*BTot/elcMass
omegaCiTot = ionCharge*BTot/ionMass
vAeTot = vAe*BTot/B0
vAiTot = vAi*BTot/B0

-- inertial length
de = vAe/omegaCe
di = vAi/omegaCi
rhoe = vtElc / omegaCe
rhoi = vtIon / omegaCi
rhoeTot = vtElc / omegaCeTot
rhoiTot = vtIon / omegaCiTot

w0 = 0.5*di
psi0 = 0.1*B0*di
guide1 = Bg
guide2 = Bg

b1 = 1.00*B0
b2 = 1.00*B0
n1 = n0/n0_ninf
n2 = n0/n0_ninf
T_e1 = elcTemp
T_i1  = ionTemp
T_i2 = T_i1

-- derived quantities (asymptotic values)
T_e2 = (0.5*(b1^2-b2^2)+0.5*(guide1^2-guide2^2)+n1*(T_i1+T_e1)-n2*T_i2)/n2 --so the system is in force balance

nuElc = 0.01*omegaCi
nuIon = nuElc/math.sqrt(ionMass)

-- domain size and simulation time
Lx = 8*math.pi*di
Ly = 4*math.pi*di
Nx = 112
Ny = 56
Ncx = 56
Ncy = 56

NvElc = 24
NvIon = 24
vMaxElc = 6*vtElc
vMaxIon = 6*vtIon
-- slightly larger velocity space extents in vz
vMaxzElcFac = 8./6.
nFieldFrames = 20
nDistFrames = 20
tFinal = 20.0/omegaCi

-- noise levels for perturbation
Bnoise_level = 0.01*B0
k0 = 1.0 --first wave mode to perturb with noise, 1.0 correspond to box size
kf = 20.0 --last wave mode to perturb with noise
Noise_index = -1.0 --spectral index of the noise

-- estimate dt and number of steps
dx = Lx/Nx
polyOrder = 2
deltaT = 2*(dx/lightSpeed)/(2*polyOrder+1)
nSteps = tFinal/deltaT

log("%50s = %g", "mi/me", ionMass / elcMass)
log("%50s = %g", "wpe/OmegaCe", wpe / omegaCe)
log("%50s = %g", "electron beta", beta_electron)
log("%50s = %g", "in-plane field", B0)
log("%50s = %g", "guide field", Bg)
log("%50s = %g", "vte/c", vtElc / lightSpeed)
log("%50s = %g", "vti/c", vtIon / lightSpeed)
log("%50s = %g", "electron plasma frequency (wpe) ", wpe)
log("%50s = %g", "electron cyclotron frequency (OmegaCe) ", omegaCe)
log("%50s = %g", "total electron cyclotron frequency (OmegaCeTot) ", omegaCeTot)
log("%50s = %g", "ion plasma frequency (wpi) ", wpi)
log("%50s = %g", "ion cyclotron frequency (OmegaCi) ", omegaCi)
log("%50s = %g", "total ion cyclotron frequency (OmegaCiTot) ", omegaCiTot)
log("%50s = %g", "electron-electron collision frequency (nuee) ", nuElc)
log("%50s = %g", "ion-ion collision frequency (nuii) ", nuIon)
log("%50s = %g", "electron inertial length (de) ", de)
log("%50s = %g", "ion inertial length (di) ", di)
log("%50s = %g", "in-plane electron gyroradius (rhoe) ", rhoe)
log("%50s = %g", "in-plane ion gyroradius (rhoi) ", rhoi)
log("%50s = %g", "total electron gyroradius (rhoeTot) ", rhoeTot)
log("%50s = %g", "total ion gyroradius (rhoiTot) ", rhoiTot)
log("%50s = %g", "Max ion velocity/vti", vMaxIon/vtIon)
log("%50s = %g", "Max elc velocity/vte", vMaxElc/vtElc)
log("%50s = %g", "Number of grid cells per di in x", Nx/(Lx/di))
log("%50s = %g", "Number of grid cells per de in x", Nx/(Lx/de))
log("%50s = %g", "tFinal ", tFinal)
log("%50s = %g", "End time in inverse ion cyclotron periods", tFinal*omegaCi)
log("%50s = %g", "Estimated time step", deltaT)
log("%50s = %g", "Estimated number of time steps", nSteps)

----------------------
--- NOISE FUNCTION ---
-----------------------
noiseGenerator = function(noiseAmp,noiseIndex,kInit,kFinal,x,y)
   local Pi = math.pi
   local sin = math.sin
   local cos = math.cos
   local sqrt = math.sqrt
   local BxNoise = 0.0
   local ByNoise = 0.0
   local JzNoise = 0.0

   local xrand_Bx = 0.0
   local xphaserand_Bx = 0.0
   local xrand_Bx2 = 0.0
   local xphaserand_Bx2 = 0.0

   local xrand_By = 0.0
   local xphaserand_By = 0.0
   local xrand_By2 = 0.0
   local xphaserand_By2 = 0.0
   local seed = 120387


   local kindex = (noiseIndex + 1.) / 2.

   for i = math.floor(kInit), math.floor(kFinal) do
       seed = seed + i
       math.randomseed(seed)
       xrand_Bx = math.random(0,2)
       xphaserand_Bx = math.random()
       xrand_Bx2 = math.random(0,2)
       xphaserand_Bx2 = math.random()

       xrand_By = math.random(0,2)
       xphaserand_By = math.random()
       xrand_By2 = math.random(0,2)
       xphaserand_By2 = math.random()


--	 For conducting y boundaries, Bx, By, and Jz = 0 at y = +/- Ly/2 and peak at y=0
--	 Az = (Lx/i*2*pi)*xrand_By*(cos(2*Pi*y/Ly)+1)^2*cos(i*2.0*Pi*x/Lx +  2*Pi*xphaserand_By)*i^kindex
	 BxNoise = BxNoise - 2.*(2.*Pi/Ly)*(Lx/(i*2.*Pi))*xrand_By*sin(2.*Pi*y/Ly)*(cos(2.*Pi*y/Ly)+1)*cos(i*2.0*Pi*x/Lx +  2.*Pi*xphaserand_By)*i^kindex
	 ByNoise = ByNoise + xrand_By*(cos(2.*Pi*y/Ly)+1)^2*sin(i*2.0*Pi*x/Lx +  2.*Pi*xphaserand_By)*i^kindex

	 JzNoise = JzNoise + (2.*Pi*i/Lx)*xrand_By*(cos(2.*Pi*y/Ly)+1)^2*cos(i*2.0*Pi*x/Lx +  2*Pi*xphaserand_By)*i^kindex + 2.*(2.*Pi/Ly)^2*(Lx/(i*2.*Pi))*xrand_By*(sin(2.*Pi*y/Ly)^2 - cos(2.*Pi*y/Ly)*(cos(2.*Pi*y/Ly)+1))*cos(i*2.0*Pi*x/Lx +  2.*Pi*xphaserand_By)*i^kindex



   end
   BxNoise = noiseAmp*BxNoise
   ByNoise = noiseAmp*ByNoise
   JzNoise = noiseAmp*JzNoise

-- This renormalizes the RMS value to noiseAmp for noiseIndex = -1
   local kdiff = math.floor(kFinal) - math.floor(kInit) + 1.0
   BxNoise = BxNoise/sqrt(2.0*kdiff*kdiff/3.0)
   ByNoise = ByNoise/sqrt(2.0*kdiff*kdiff/3.0)
   JzNoise = JzNoise/sqrt(2.0*kdiff*kdiff/3.0)

   return BxNoise, ByNoise, JzNoise
end



-- Maxwellian in 2x3v
local function maxwellian3D(n, vx, vy, vz, ux, uy, uz, temp, mass)
   local v2 = (vx - ux)^2 + (vy - uy)^2 + (vz - uz)^2
   return n/math.sqrt((2*math.pi*temp/mass)^3)*math.exp(-mass*v2/(2.0*temp))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = tFinal, -- end time
   nFrame = nFieldFrames, -- number of output frames
   lower = {-Lx/2., -Ly/2.}, -- configuration space lower left
   upper = {Lx/2., Ly/2.}, -- configuration space upper right
   cells = {Nx, Ny}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = polyOrder, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   parallelizeSpecies = true,
   -- decomposition for configuration space
   decompCuts = {Ncx, Ncy}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory
   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions
   
   restartFrameEvery = 0.2,
   -- integrated moment flag, compute quantities 1000 times in simulation
   calcIntQuantEvery = 0.001,

   -- electrons
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-vMaxElc, -vMaxElc, -vMaxElc*vMaxzElcFac},
      upper = {vMaxElc, vMaxElc, vMaxElc*vMaxzElcFac},
      cells = {NvElc, NvElc, NvElc},
      -- initial conditions
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local tanh = math.tanh
         local cosh = math.cosh
         local cos = math.cos
         local sin = math.sin

         local me = elcMass
         local mi = ionMass
         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

	 
         local b1x = -0.5*(b2+b1)*(tanh(y/w0)-tanh((y+Ly)/w0) - tanh((y-Ly)/w0))
         local b1z = Bg
         local b1y = 0.0
	 
  	 local Tvali = 0.5*(T_i2+T_i1)
         local Tvale = 0.5*(T_e2+T_e1)
         local n = (0.5*(b1^2 - b1x^2) + n1*(T_i1+T_e1))/(Tvali+Tvale) -- only correct for mu0 = 1


         local TeFrac = Tvale/(Tvale + Tvali)
         local TiFrac = Tvali/(Tvale + Tvali)

         local Jx = 0.0
         local Jz  = 0.5*(b2+b1)/w0*((1.0/cosh(y/w0))^2 - (1.0/cosh((y+Ly)/w0))^2 - (1.0/cosh((y-Ly)/w0))^2) - psi0*cos(_2pi*x/Lx)*cos(Pi*y/Ly)*((_2pi/Lx)^2 + (Pi/Ly)^2)
	 local BxNoise, ByNoise, JzNoise = noiseGenerator(Bnoise_level, Noise_index, k0, kf, x, y)	 
	 Jz = Jz + JzNoise

         local Jxe = Jx*TeFrac
         local Jxi = Jx*TiFrac
         local Jze = Jz*TeFrac
         local Jzi = Jz*TiFrac

         local vdrift_x = Jxe/(elcCharge*n)
         local vdrift_y = 0.0
         local vdrift_z = Jze/(elcCharge*n)

	 local fv = maxwellian3D(n, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, elcTemp, elcMass)
       
	 return fv
      end,
     -- boundary conditions
      bcy = { Plasma.ReflectBC{}, Plasma.ReflectBC{} },
      nDistFuncFrame = nDistFrames,
      evolve = true, -- evolve species?
      diagnostics = { "M0", "M1i", "M2", "M2ij", "M3i", "intM0", "intM1i", "intM2Flow", "intM2Thermal" }, 
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },
   -- protons
   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-vMaxIon, -vMaxIon, -vMaxIon},
      upper = {vMaxIon, vMaxIon, vMaxIon},
      cells = {NvIon, NvIon, NvIon},
      -- initial conditions
      init = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]

         local tanh = math.tanh
         local cosh = math.cosh
         local cos = math.cos
         local sin = math.sin


         local me = elcMass
         local mi = ionMass
         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi


         local b1x = -0.5*(b2+b1)*(tanh(y/w0)-tanh((y+Ly)/w0) - tanh((y-Ly)/w0))
         local b1z = Bg
         local b1y = 0.0
	 
  	 local Tvali = 0.5*(T_i2+T_i1)
         local Tvale = 0.5*(T_e2+T_e1)
         local n = (0.5*(b1^2 - b1x^2) + n1*(T_i1+T_e1))/(Tvali+Tvale) -- only correct for mu0 = 1


         local TeFrac = Tvale/(Tvale + Tvali)
         local TiFrac = Tvali/(Tvale + Tvali)

         local Jx = 0.0
         local Jz  = 0.5*(b2+b1)/w0*((1.0/cosh(y/w0))^2 - (1.0/cosh((y+Ly)/w0))^2 - (1.0/cosh((y-Ly)/w0))^2) - psi0*cos(_2pi*x/Lx)*cos(Pi*y/Ly)*((_2pi/Lx)^2 + (Pi/Ly)^2)
	 local BxNoise, ByNoise, JzNoise = noiseGenerator(Bnoise_level, Noise_index, k0, kf, x, y)	 
	 Jz = Jz + JzNoise

         local Jxe = Jx*TeFrac
         local Jxi = Jx*TiFrac
         local Jze = Jz*TeFrac
         local Jzi = Jz*TiFrac

         local vdrift_x = Jxi/(ionCharge*n)
         local vdrift_y = 0.0
         local vdrift_z = Jzi/(ionCharge*n)


	 local fv = maxwellian3D(n, vx, vy, vz, vdrift_x, vdrift_y, vdrift_z, ionTemp, ionMass)

	 return fv
      end,
      -- boundary conditions
      bcy = { Plasma.ReflectBC{}, Plasma.ReflectBC{} },
      nDistFuncFrame = nDistFrams,
      evolve = true, -- evolve species?
      diagnostics = { "M0", "M1i", "M2", "M2ij", "M3i", "intM0", "intM1i", "intM2Flow", "intM2Thermal" }, 
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
   },

   -- field solver
   field = Plasma.Field {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         
         local tanh = math.tanh
         local cosh = math.cosh
         local cos = math.cos
         local sin = math.sin

         local me = elcMass
         local mi = ionMass
         local qe = elcCharge
         local qi = ionCharge

         local Pi = math.pi
         local _2pi = 2.0*Pi
         local _4pi = 2.0*_2pi

	 
         local b1x = -0.5*(b2+b1)*(tanh(y/w0)-tanh((y+Ly)/w0) - tanh((y-Ly)/w0))
         local b1z = Bg
         local b1y = 0.0
	 
     	 local Tvali = 0.5*(T_i2+T_i1)
         local Tvale = 0.5*(T_e2+T_e1)
         local n = (0.5*(b1^2 - b1x^2) + n1*(T_i1+T_e1))/(Tvali+Tvale)

	 local BxNoise, ByNoise, JzNoise = noiseGenerator(Bnoise_level, Noise_index, k0, kf, x, y)	 
	 local Bx = b1x + psi0*Pi/Ly*cos(_2pi*x/Lx)*sin(Pi*y/Ly) + BxNoise
         local By = b1y - psi0*_2pi/Lx*sin(_2pi*x/Lx)*cos(Pi*y/Ly) + ByNoise
         local Bz = b1z

         local TeFrac = Tvale/(Tvale + Tvali)
         local TiFrac = Tvali/(Tvale + Tvali)

         local Jx = 0.0
	 local Jz  = 0.5*(b2+b1)/w0*((1.0/cosh(y/w0))^2 - (1.0/cosh((y+Ly)/w0))^2 - (1.0/cosh((y-Ly)/w0))^2) - psi0*cos(_2pi*x/Lx)*cos(Pi*y/Ly)*((_2pi/Lx)^2 + (Pi/Ly)^2)
         Jz = Jz + JzNoise

         local Jxe = Jx*TeFrac
         local Jxi = Jx*TiFrac
         local Jze = Jz*TeFrac
         local Jzi = Jz*TiFrac

         -- Assumes qi = abs(qe)
         local u_xe = Jxe/(elcCharge*n)
         local u_ye = 0.0
         local u_ze = Jze/(elcCharge*n)
   
         -- E = - v_e x B ~  (J - u) x B
         local Ex = - 0.0*(u_ye*Bz - u_ze*By)
         local Ey = - 0.0*(u_ze*Bx - u_xe*Bz)
         local Ez = - 0.0*(u_xe*By - u_ye*Bx)

	 return Ex, Ey, Ez, Bx, By, Bz
      end,
      bcy = { Plasma.Field.bcReflect, Plasma.Field.bcReflect },
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
