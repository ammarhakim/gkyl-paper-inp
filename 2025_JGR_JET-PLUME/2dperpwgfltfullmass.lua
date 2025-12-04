-- Gkeyll
-- Continuum kinetic 2X2V simulation of the weibel instability with guiding field
-- Collin Brown

local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()

----------------------------------------------------------------------
pi = math.pi

-- Constants
chargeElc = -1.0
chargeIon = 1.0
massElc = 1.0
massIon = 1836.0
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1.0/math.sqrt(epsilon0*mu0)

-- Initial conditions
nElc10 = 0.5
nElc20 = 0.5
uxElc10 = 0.0
uyElc10 = 0.3
uzElc10 = 0.0
uxElc20 = 0.0
uyElc20 = -0.3
uzElc20 = 0.0
TElc10 = 0.001
TElc20 = 0.001
nIon0 = 1.0
uxIon0 = 0.0
uyIon0 = 0.0
uzIon0 = 0.0
TIon0 = 0.001
k0 = 0.4
perturb = 1e-3

vthElc10 = math.sqrt(TElc10/massElc)
vthElc20 = math.sqrt(TElc20/massElc)
vthIon0 = math.sqrt(TIon0/massIon)

polyOrder = 2
numCellsConf = {32,32}
numCellsPhase = {64, 64} --Nvx, Nvy
lowerBndConf = {0.0,0.0}
upperBndConf = {6*pi/k0,6*pi/k0}
lowerBndPhaseElc = {-vthElc10*50, -vthElc10*50} --vx, vy
upperBndPhaseElc = {vthElc10*50, vthElc10*50} --vx, vy
lowerBndPhaseIon = {-20*vthIon0, -20*vthIon0} --vx, vy
upperBndPhaseIon = {20*vthIon0, 20*vthIon0} --vx, vy
periodicDirs = {1,2}
cfl = (1.0/3.0)/(2*polyOrder + 1)
tStart = 0.0
tEnd = 300.
numFrames = 200

function maxwellian2D(n, vx, vy, ux, uy, vth)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2*pi*vth^2)*math.exp(-v2/(2*vth^2))
end

vlasovApp = Vlasov.App {
  --------------------------------------------------------------------------------
  -- Common
  --------------------------------------------------------------------------------
  logToFile = true,

  tEnd = tEnd,                   -- End time
  nFrame = numFrames,            -- Number of output frames
  lower = lowerBndConf,          -- Lower boundary of configuration space
  upper = upperBndConf,          -- Upper boundary of configuration space
  cells = numCellsConf,          -- Configuration space cells
  basis = "serendipity",         -- One of "serendipity", "maximal-order", or "tensor"
  polyOrder = polyOrder,         -- Polynomial order
  timeStepper = "rk3",           -- One of "rk2", "rk3", or "rk3s4"
  periodicDirs = periodicDirs,
  decompCuts = {16,1}, 
  
  elc = Vlasov.Species {
     charge = chargeElc, mass = massElc,
     -- Velocity space grid
     lower = {lowerBndPhaseElc[1],lowerBndPhaseElc[2]},
     upper = {upperBndPhaseElc[1],upperBndPhaseElc[2]},
     cells = {numCellsPhase[1], numCellsPhase[2]},
     -- Initial conditions
     init = function (t, xn)
        local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
        
        return maxwellian2D(nElc10, vx, vy, uxElc10, uyElc10, vthElc10) +  maxwellian2D(nElc20, vx, vy, uxElc20, uyElc20, vthElc20)

     end,
     evolve = true,
  },
  
  -- ions
  ions = Vlasov.Species {
     charge = chargeIon, mass = massIon,
     -- velocity space grid
     lower = {lowerBndPhaseIon[1],lowerBndPhaseIon[2]},
     upper = {upperBndPhaseIon[1],upperBndPhaseIon[2]},
     cells = {numCellsPhase[1], numCellsPhase[2]},

     -- initial conditions
     init = function (t, xn)
        local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
       
        return maxwellian2D(nIon0, vx, vy, uxIon0, uyIon0, vthIon0)
     end,
     evolve = true, -- evolve species?
  },  
  
  --------------------------------------------------------------------------------
  -- Field solver
  --------------------------------------------------------------------------------
  field = Vlasov.Field {
     epsilon0 = epsilon0, mu0 = mu0,
     -- hasMagneticField = false,
     init = function (t, xn)
        local x,y = xn[1],xn[2]
        local Ex = 0.0
        local Ey = perturb*math.sin(k0*(x)) 
        local Ez = 0.0
        local Bx = 0.0  --keep guiding field in plane of simulation
        local By = 0.2  --keep guiding field in plane of simulation
        local Bz = 0.0
        return Ex, Ey, Ez, Bx, By, Bz
     end,
     evolve = true,
   },

}

--------------------------------------------------------------------------------
-- Run application
--------------------------------------------------------------------------------
vlasovApp:run()
  
  
