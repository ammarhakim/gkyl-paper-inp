-- Gkyl ------------------------------------------------------------------------
-- Simulation based off data from Chen and Wang at AGU 2018 --------------------
-- Nearly transverse shock, approximated as exactly transverse -----------------
-- shock velocity lowered to reduce computational cost -------------------------
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell
local Constants = require "Lib.Constants"
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = True
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

-- global parameters
gasGamma = 5.0/3.0
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1.0

--Assume MMS Plasma is Hydrogen 
elcCharge = -1.0
ionCharge = 1.0
elcMass = 1.0
ionMass = 100.0

--plasma parameters
n0 = 1.0
beta_electron = 0.7
beta_proton = 1.3
vte = 1.0/16.0 --vte/c = 1/8
Ti_Te = beta_proton/beta_electron
vti = vte*math.sqrt(Ti_Te)/math.sqrt(ionMass/elcMass)

Te = vte^2*elcMass/2.0
Ti = vti^2*ionMass/2.0

vAe = vte/math.sqrt(beta_electron)
vA = vti/math.sqrt(beta_proton)

--Upstream magnetic field is vAe in normalized units
B0 = vAe

--Derived parameters
cs = math.sqrt(gasGamma*(Te+Ti)/ionMass)

wpe = math.sqrt(ionCharge^2*n0/(epsilon0*elcMass))
wpi = math.sqrt(ionCharge^2*n0/(epsilon0*ionMass))

OmegaCi = ionCharge*B0/ionMass
OmegaCe = ionCharge*B0/elcMass

di = lightSpeed/wpi
de = lightSpeed/wpe

rhoi = vti/OmegaCi
rhoe = vte/OmegaCe

lambdaD = vte/wpe/math.sqrt(2)

--initial shock speed
u_shock = 3.0*vA
--small noise to Bz to excite parallel/transverse dynamics
pert = 1.0e-3

nuElc = 0.01*OmegaCi
nuIon = nuElc/math.sqrt(ionMass)

polyOrder = 2
xlower = 0.0
xupper = 25*di
nx = 1536

--velocity space extents
vMaxElc = 8.0*vte
vMaxIon = 16.0*vti
nvElc = 64
nvIon = 64

-- determine how long to run the simulation
tEnd = 11.0/OmegaCi
nFrames = 110

-- estimate dt and number of steps
dx = (xupper-xlower)/nx
deltaT = 2*(dx/lightSpeed)/(2*polyOrder+1)
nSteps = tEnd/deltaT

log("%50s = %g", "mi/me", ionMass / elcMass)
log("%50s = %g", "wpe/OmegaCe", wpe / OmegaCe)
log("%50s = %g", "electron beta", beta_electron)
log("%50s = %g", "proton beta", beta_proton)
log("%50s = %g", "vte/c", vte / lightSpeed)
log("%50s = %g", "vti/c", vti / lightSpeed)
log("%50s = %g", "electron plasma frequency (wpe) ", wpe)
log("%50s = %g", "electron cyclotron frequency (OmegaCe) ", OmegaCe)
log("%50s = %g", "ion plasma frequency (wpi) ", wpi)
log("%50s = %g", "ion cyclotron frequency (OmegaCi) ", OmegaCi)
log("%50s = %g", "electron inertial length (de) ", de)
log("%50s = %g", "ion inertial length (di) ", di)
log("%50s = %g", "electron Debye length (lambdaD) ", lambdaD)
log("%50s = %g", "time of shock propagation ", xupper/u_shock)
log("%50s = %g", "Sonic Mach number", u_shock/cs)
log("%50s = %g", "Alfven Mach number", u_shock/vA)
log("%50s = %g", "Shock velocity/vti", u_shock/vti)
log("%50s = %g", "Shock velocity/vte", u_shock/vte)
log("%50s = %g", "Max ion velocity/vti", vMaxIon/vti)
log("%50s = %g", "Max elc velocity/vte", vMaxElc/vte)
log("%50s = %g", "Number of grid cells per di in x", nx/((xupper-xlower)/di))
log("%50s = %g", "Number of grid cells per de in x", nx/((xupper-xlower)/de))
log("%50s = %g", "tEnd ", tEnd)
log("%50s = %g", "End time in inverse ion cyclotron periods", tEnd*OmegaCi)
log("%50s = %g", "Estimated time step", deltaT)
log("%50s = %g", "Estimated number of time steps", nSteps)

-- initialization function
local function maxwellian2D(n, vx, ux, vy, uy, temp, mass)
   v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*math.pi*temp/mass)*math.exp(-mass*v2/(2*temp))
end

sim = Plasma.App {
   logToFile = true,

   tEnd = tEnd, -- end time
   nFrame = nFrames, -- number of output frames
   lower = {xlower}, -- configuration space lower left
   upper = {xupper}, -- configuration space upper right
   cells = {nx}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = polyOrder, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {32}, -- cuts in each configuration direction
   useShared = true, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Plasma.Species {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-vMaxElc, -vMaxElc},
      upper = {vMaxElc, vMaxElc},
      cells = {nvElc, nvElc},
      -- initial conditions
      init = function (t, xn)
         local x, vx, vy = xn[1], xn[2], xn[3]
         local fv = maxwellian2D(n0, vx, -u_shock, vy, 0.0, Te, elcMass)
         return fv
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcReflect, Plasma.Species.bcCopy },
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      coll = Plasma.LBOCollisions {
         collideWith = {'elc'},
         frequencies = {nuElc},
      },
   },

   ion = Plasma.Species {
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-vMaxIon, -vMaxIon},
      upper = {vMaxIon, vMaxIon},
      cells = {nvIon, nvIon},
      -- initial conditions
      init = function (t, xn)
         local x, vx, vy = xn[1], xn[2], xn[3]
         local fv = maxwellian2D(n0, vx, -u_shock, vy, 0.0, Ti, ionMass)
         return fv
      end,
      evolve = true, -- evolve species?
      bcx = { Plasma.Species.bcReflect, Plasma.Species.bcCopy },
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      coll = Plasma.LBOCollisions {
         collideWith = {'ion'},
         frequencies = {nuIon},
      },
   },
   
   -- field solver
   field = Plasma.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x = xn[1]
         local Ex = 0.0
         local Ey = -u_shock*B0
         local Ez = 0.0
         local Bx = 0.0
         local By = 0.0
         local Bz = B0
         return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- evolve field?
      bcx = { Plasma.Field.bcReflect, Plasma.Field.bcCopy },
   },
}
-- run application
sim:run()
