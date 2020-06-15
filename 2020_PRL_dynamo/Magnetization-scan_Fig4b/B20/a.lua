-- Gkyl ------------------------------------------------------------------------

local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell
local Constants = require "Lib.Constants" --contains universal constants taken from NIST's website
local Logger = require "Lib.Logger"
local prng = require "sci.prng"

local logger = Logger {
   logToFile = true
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

-- ******************************
--     INPUT PARAMETERS
-- ******************************


Cnu = 0.05
BoxMultiplier = 8.0 

-- normalization parameters, working in SI units

epsilon0 =  Constants.EPSILON0 -- permittivity of free space

mu0 = Constants.MU0 -- pemiability of free space

lightSpeed = Constants.SPEED_OF_LIGHT -- speed of light

unitCharge = Constants.ELEMENTARY_CHARGE

elcMass = Constants.ELECTRON_MASS -- electron mass

elcCharge = -1.0 * unitCharge -- electron charge

ionMass = Constants.PROTON_MASS -- assuming proton mass

ionCharge = 1.0 * unitCharge -- assuming proton charge

eVJ = 1.0 * unitCharge -- 1 eV in Joules 

-- input quantities

TeeV = 1000  -- Te in eV

TieV = 1000  -- Ti in eV

M0 = 0.5  --  Mach number at scale L0. u0 = M*cs

ne = 2.302e+28 -- electron density

k0 = 822500.0/BoxMultiplier -- inverse of the length-scale L0 

lnLambda = 10 -- Coulomb logarithm 

-- Noise generation parameters

BNoiseAmpl = 20.0 -- Amplitude of B noise (in Tesla)
 
nmodes = 1 -- Number of modes (mode 1 has a wave length comparable to the box)

-- Forcing related quantities (Here only ions are forced)

forcefactor = 3.0 -- Forcing is set to be Force = forceFactor * u * m / tstreamIon, 
                  -- (note that due to normalization m is not in the forcing function)
                  -- We assume that the ions loose their momentum on the ion streaming time scale  

-- derived quantities

Te = eVJ*TeeV -- Te in Joules 

Ti = eVJ*TieV -- Ti in Joules

vte = math.sqrt(2*Te / elcMass) -- electron thermal speed 

vti = math.sqrt(2*Ti / ionMass) -- ion thermal speed

wpe = math.sqrt(ionCharge^2*ne/(epsilon0*elcMass)) -- electron plasma frequency

wpi = math.sqrt(ionCharge^2*ne/(epsilon0*ionMass)) -- ion plasma frequency

OmegaCi = ionCharge*BNoiseAmpl/ionMass -- ion cyclotron frequency from B noise

OmegaCe = ionCharge*BNoiseAmpl/elcMass -- electron cyclotron frequency from B noise

di = lightSpeed/wpi -- ion inertial length

de = lightSpeed/wpe -- electron inertial length

lambdaD = vte/wpe/math.sqrt(2)

cs = math.sqrt(Te / ionMass) -- proxy for sound speed

u0 = cs*M0 -- flow velocity size at scale L0 

tauei = 47.2488 * math.sqrt(elcMass)*(Te)^(3/2)*epsilon0^2/(ne*unitCharge^4*lnLambda) -- e-i collision time  

nuei = Cnu* 1.0/tauei -- collision frequency is inverse of collision time

nuee = Cnu* 2.0*nuei -- estimating electron-electron collision time from electron-ion collision time

sigma = 1.96928 * ne * unitCharge^2 * tauei / elcMass -- Spitzer conductivity

eta = 1/(sigma*mu0) -- Magnetic diffusivity

nuii = Cnu* ne*unitCharge^4*lnLambda/(4*math.pi*epsilon0^2*ionMass^2*vti^3) -- i-i collision frequency

nuie = Cnu* nuei*elcMass/ionMass -- estimating ion-electron collision frequency from ion-electron collision frequency

nu = 1.80477*Ti/(nuii*ionMass) -- kinematic viscosity

Re = u0 / (nu*k0) -- Reynolds number at scale L0

Rm = u0 / (eta*k0) -- Magnetic Reynolds number at scale L0

L0 = 1/k0 -- System size

tstreamIon = L0/vti -- Thermal ion streaming time across simulation domain 

tTurnover = L0/u0 -- Eddy turnover time on the largest (box-size) scale

tLightCrossing = L0/lightSpeed -- Light crossing time for simulation domain 

tEnd = 1.0 * tstreamIon -- end of simulation time (now: fraction of the ion thermal streaming time)

rLe = - vte / (elcCharge * BNoiseAmpl / elcMass) -- thermal electron Larmor radius at B noise amplitude

rLi = - vti / (ionCharge * BNoiseAmpl / ionMass) -- thermal ion Larmor radius at B noise amplitude

-- domain size and simulation time

LX, LY, LZ = 1*L0, 1*L0, 1*L0

NX, NY, NZ = 40, 40, 40

NVX, NVY, NVZ = 20, 20, 20

viLimits, veLimits = 3.0, 3.0

-- estimate dt and number of steps
dx = LX/NX
deltaT = 0.9*(dx/lightSpeed)
nSteps = tEnd/deltaT

-- Printing some information

log("%50s = %g", "mi/me", ionMass / elcMass)
log("%50s = %g", "wpe/OmegaCe", wpe / OmegaCe)
log("%50s = %g", "Electron temperature in Joule", Te)
log("%50s = %g", "ion temperature in Joule", Ti)
log("%50s = %g", "Electron temperature in electron-volts", TeeV)
log("%50s = %g", "ion temperature in electron-volts", TieV)
log("%50s = %g", "vte/c", vte / lightSpeed)
log("%50s = %g", "vti/c", vti / lightSpeed)
log("%50s = %g", "electron plasma frequency (wpe) in Hz", wpe)
log("%50s = %g", "electron cyclotron frequency (OmegaCe) in Hz", OmegaCe)
log("%50s = %g", "ion plasma frequency (wpi) in Hz", wpi)
log("%50s = %g", "ion cyclotron frequency (OmegaCi) in Hz", OmegaCi)
log("%50s = %g", "electron inertial length (de) in meters", de)
log("%50s = %g", "ion inertial length (di) in meters", di)
log("%50s = %g", "electron Debye length (lambdaD) in meters", lambdaD)
log("%50s = %g", "Mach number", u0/cs)
log("%50s = %g", "Shock velocity/vti", u0/vti)
log("%50s = %g", "Shock velocity/vte", u0/vte)
log("%50s = %g", "Number of grid cells per di in x", NX/(LX/di))
log("%50s = %g", "Number of grid cells per de in x", NX/(LX/de))
log("%50s = %g", "Number of grid cells per lambdaD in x", NX/(LX/lambdaD))
log("%50s = %g", "Number of grid cells per rho_e in x", NX/(LX/rLe))
log("%50s = %g", "rho_e per domain size", rLe/LX)
log("%50s = %g", "tEnd in seconds", tEnd)
log("%50s = %g", "Estimated time step", deltaT)
log("%50s = %g", "Estimated number of time steps", nSteps)

-- ******************************
--       FUNCTIONS
-- ******************************

-- Maxwellian in 3 velocity dimensions (assumes temperature in Joules)

local function maxwellian3D(n, vx, vy, vz, ux, uy, uz, mass, temp)

   local v2 = (vx - ux)^2 + (vy - uy)^2 + (vz - uz)^2 

   return n*(mass/(2*math.pi*temp))^(3/2)*math.exp(-mass*v2/(2*temp))

end

local function maxwellian2D(n, vx, vy, ux, uy, mass, temp)

   local v2 = (vx - ux)^2 + (vy - uy)^2 

   return n*(mass/(2*math.pi*temp))*math.exp(-mass*v2/(2*temp))

end

-- Noise generator

local function noiseGenerator(BNoiseAmpl,nmodes,x,y,z)

   local Pi = math.pi
   local _2pi = 2.0 * math.pi
   local sin = math.sin
   local cos = math.cos
   local sqrt = math.sqrt

   local Bx = 0.0
   local By = 0.0
   local Bz = 0.0

   local Jx = 0.0
   local Jy = 0.0
   local Jz = 0.0

   local B_xy = 0.0
   local B_xz = 0.0
   local B_yx = 0.0
   local B_yz = 0.0
   local B_zx = 0.0
   local B_zz = 0.0

   local phase_xy = 0.0
   local phase_xz = 0.0
   local phase_yx = 0.0
   local phase_yz = 0.0
   local phase_zx = 0.0
   local phase_zz = 0.0

   local seed = 120387

   for i = 1, nmodes do
      seed = seed + 1
      math.randomseed(seed)
      B_xy = math.random()
      B_xz = math.random()
      B_yx = math.random()
      B_yz = math.random()
      B_zx = math.random()
      B_zy = math.random()
      phase_xy = math.random()
      phase_xz = math.random()
      phase_yx = math.random()
      phase_yz = math.random()
      phase_zx = math.random()
      phase_zy = math.random()

      Bx = Bx + B_xy*cos(_2pi*i*(y/LY+phase_xy)) + B_xz*cos(_2pi*i*(z/LZ+phase_xz))
      By = By + B_yx*cos(_2pi*i*(x/LX+phase_yx)) + B_yz*cos(_2pi*i*(z/LZ+phase_yz))
      Bz = Bz + B_zx*cos(_2pi*i*(x/LX+phase_zx)) + B_zy*cos(_2pi*i*(y/LY+phase_zy))

      Jx = Jx + i*( B_yz*sin( _2pi*i*(z/LZ+phase_yz) )/LZ - B_zy*sin( _2pi*i*(y/LY+phase_zy) )/LY )
      Jy = Jy - i*( B_xz*sin( _2pi*i*(z/LZ+phase_xz) )/LZ + B_zx*sin( _2pi*i*(x/LX+phase_zx) )/LX )
      Jz = Jz + i*( B_xy*sin( _2pi*i*(y/LY+phase_xy) )/LY - B_yx*sin( _2pi*i*(x/LX+phase_yx) )/LX )

   end
      Bx = Bx * BNoiseAmpl
      By = By * BNoiseAmpl
      Bz = Bz * BNoiseAmpl
      Jx = Jx * BNoiseAmpl * _2pi/mu0
      Jy = Jy * BNoiseAmpl * _2pi/mu0
      Jz = Jz * BNoiseAmpl * _2pi/mu0

   return Bx, By, Bz, Jx, Jy, Jz
end

local function BJGenerator(BNoiseAmpl,nmodes,x)

   local Pi = math.pi
   local _2pi = 2.0 * math.pi
   local sin = math.sin
   local cos = math.cos
   local sqrt = math.sqrt

   local Bx = 0.0
   local By = 0.0
   local Bz = 0.0

   local Jx = 0.0
   local Jy = 0.0
   local Jz = 0.0

   local B_xy = 0.0
   local B_xz = 0.0
   local B_yx = 0.0
   local B_yz = 0.0
   local B_zx = 0.0
   local B_zz = 0.0

   local phase_xy = 0.0
   local phase_xz = 0.0
   local phase_yx = 0.0
   local phase_yz = 0.0
   local phase_zx = 0.0
   local phase_zz = 0.0

   local seed = 120387

   for i = 1, nmodes do
      seed = seed + 1
      math.randomseed(seed)
      B_zx = 1
      phase_zx = math.random()

      Bz = Bz + B_zx*cos(_2pi*i*(x/LX+phase_zx))
      Jy = Jy + i*B_zx*sin( _2pi*i*(x/LX+phase_zx) )/LX

   end
      Bz = Bz * BNoiseAmpl
      Jy = Jy * BNoiseAmpl * _2pi/mu0

   return Bx, By, Bz, Jx, Jy, Jz
end


-- Generates Roberts flow 1 on cubic domain, flow speed is called u0 

local function robertsFlow(x,y,z)  -- Assumes that LX=LY=LZ and the flow speed u0 are defined

   local Pi = math.pi
   local _2pi = 2.0 * math.pi
   local sin = math.sin
   local cos = math.cos   

   local uRx, uRy, uRz = 0.0, 0.0, 0.0

   uRx = u0 * ( cos(_2pi*y/LY) - cos(_2pi*z/LZ) )
   uRy = u0 * sin(_2pi*z/LZ)
   uRz = u0 * sin(_2pi*y/LY)
   return uRx, uRy, uRz
end  

-- ******************************
--   VLASOV APP
-- ******************************

vlasovApp = Plasma.App {
   logToFile = true,

   tEnd = tEnd, -- end time 

   suggestedDt = deltaT, -- suggested time step (here: some fraction of the light crossing time across cell)

   nFrame = 20, -- number of output frames

   lower = {0.0}, -- configuration space lower left

   upper = {LX}, -- configuration space upper right

   cells = {NX}, -- configuration space cells

   basis = "serendipity", -- one of "serendipity" or "maximal-order"

   polyOrder = 2, -- polynomial order

   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- restartFrameEvery = 0.1, -- not using restart frames for the dime being

   -- decomposition for configuration space

   decompCuts = {1}, -- cuts in each configuration direction

   useShared = true, -- if to use shared memory

   periodicDirs ={1},

   -- electrons

   elc = Plasma.Species {
      
      nDistFuncFrame = 10,

      charge = elcCharge, mass = elcMass,

      -- velocity space grid

      lower = {-veLimits*vte, -veLimits*vte},

      upper = {veLimits*vte, veLimits*vte},

      cells = {NVX, NVY},

      decompCuts = {1, 1}, -- do not change, no parallelization in velocity space currently

      -- initial conditions

      init = function (t, xn)

	 local x, vx, vy= xn[1], xn[2], xn[3]
         
         local ux, uy = 0.0, 0.0

         local Bx, By, Bz, Jx, Jy, Jz = BJGenerator(BNoiseAmpl,nmodes,x)

         uy = Jy/(elcCharge*ne)   

	 local fv = maxwellian2D(ne, vx, vy, ux, uy, elcMass, Te)

	 return fv

      end,

      evolve = true, -- evolve species?

      -- write out density, flow, total energy, and heat flux moments

      diagnosticMoments = { "M0", "M1i", "M2", "M2ij" },

      -- Collisions.
      coll = Plasma.LBOCollisions {

	 collideWith  = { "elc", "ion"},

     	 frequencies  = { nuee, nuei},

         -- Optional arguments:

         -- crossOption = "Greene",    -- Or crossOption="HeavyIons".

         -- betaGreene  = 1.0,
      },
   },

   -- protons

   ion = Plasma.Species {
            
      nDistFuncFrame = 10,

      charge = ionCharge, mass = ionMass,

      -- velocity space grid

      lower = {-viLimits*vti, -viLimits*vti},

      upper = {viLimits*vti, viLimits*vti},

      cells = {NVX, NVY},

      decompCuts = {1, 1}, -- do not change, no parallelization in velocity space currently

      -- initial conditions

      init = function (t, xn)

	 local x, vx, vy = xn[1], xn[2], xn[3]
         
         local ux, uy = 0.0, 0.0

	 local fv = maxwellian2D(ne, vx, vy, ux, uy, ionMass, Ti)

	 return fv

      end,

      -- Forcing (set to be proportional to the Roberts flow.)
      -- vlasovExtForceFunc = function(t, xn)

         -- x, y, z = xn[1], xn[2], xn[3]
         
         -- local ux, uy, uz = 0.0, 0.0, 0.0

         -- local uRx, uRy, uRz = robertsFlow(x,y,z)
         
         -- force_x = forcefactor*uRx / tstreamIon 
 
         -- force_y = forcefactor*uRy / tstreamIon    

         -- force_z = forcefactor*uRz / tstreamIon      
         
         -- return force_x, force_y, force_z
      -- end,

      evolve = true, -- evolve species?

      -- write out density, flow, total energy, and heat flux moments

      diagnosticMoments = { "M0", "M1i", "M2", "M2ij"},

      -- Collisions.
      coll = Plasma.LBOCollisions {

	 collideWith = { "ion", "elc"},

     	 frequencies = { nuii, nuie},

         -- Optional arguments:

         --crossOption = "Greene",    -- Or crossOption="HeavyIons".

         --betaGreene  = 1.0,
      },
   },

   -- field solver

   field = Plasma.Field {

      epsilon0 = epsilon0, mu0 = mu0,

      init = function (t, xn)

	 local x, y, z = xn[1], xn[2], xn[3]

         local Bx, By, Bz, Jx, Jy, Jz = BJGenerator(BNoiseAmpl,nmodes,x)

         local Ex, Ey, Ez = 0.0, 0.0, 0.0

	 return Ex, Ey, Ez, Bx, By, Bz

      end,

      evolve = true, -- evolve field?

   },

-- ******************************
--   EXECUTING VLASOV APP
-- ******************************

}
-- run application
vlasovApp:run()
