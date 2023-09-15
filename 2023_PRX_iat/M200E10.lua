-- Gkyl ------------------------------------------------------------------------
-- Z.Liu 6/25/2022
-- Mass=25, Temp=50
-- E1, E_ext = 1.0e-5
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- electron parameters
vDriftElc = {0.0,0.0}   ---modified from 0.159
vtElc = 0.02
-- ion parameters
vDriftIon = {0.0,0.0}
vtIon = 0.0002      ---- modified from 0.001 (use Te = 50 Ti)
tempRatio = 50.0
-- mass ratio
massRatio = 200.0  ----modified from 25 

perturbation = 1.0e-4 -- distribution function perturbation
noise = 1.0e-6

omegaPe = 1.0
omegaPi = omegaPe/math.sqrt(massRatio)
lambdaDe = vtElc / omegaPe

-- Collision frequencies.
nuee = 0.0001
nuei = 0.0 --nuee              --RLW 
nuii = 2.0e-5  --should have been nuee/math.sqrt(massRatio)*math.pow(tempRatio,1.5)
nuie = 0.0

lx = {1.0,0.5}
nx = {32,16}

pOrder = 2

accion = -2.5e-7

local function maxwellian2v(v, vDrift, vt)
    if v[1]-vDrift[1] > 4.0*vtElc then	
       return 0.00001*1/(2*math.pi*vt^2)*math.exp(-((v[1]-vDrift[1])^2+(v[2]-vDrift[2])^2)/(2*vt^2))
    elseif v[1]-vDrift[1] <- 4.0*vtElc then
       return 0.00001*1/(2*math.pi*vt^2)*math.exp(-((v[1]-vDrift[1])^2+(v[2]-vDrift[2])^2)/(2*vt^2))
    else
       return 1/(2*math.pi*vt^2)*math.exp(-((v[1]-vDrift[1])^2+(v[2]-vDrift[2])^2)/(2*vt^2))
    end
end


local function sponEmissionSource(x_table, t, lx_table, ncells_table, p)
    -- x_table = {x, y} are independent variables.
    -- t is the time
    -- lx_table = {xmax-xmin, ymax - ymin}
    -- ncells_table = {ncellsx, ncellsy}
    -- p is polynomial order.
    ------------------------
    -- returns stochastic source term (\delta E_x, \delta E_y, delta f_i, \delta f_e) that models spontaneous emission of ion-acoustic waves
    local Nx = ncells_table[1]*(p+1) -- Number of spatial degrees of freedom along x
    local Ny = ncells_table[2]*(p+1) -- Number of spatial degrees of freedom along y
    local Lx, Ly = lx_table[1], lx_table[2]
    local x, y, vx, vy = x_table[1], x_table[2], x_table[3], x_table[4]
    local fIon = 0.0
    local ksquared = 0.0
    math.randomseed(math.floor(1000000*t)) --I want all xs and ys (and vx, vy) to see the same random phase at a given time step.  Since t will be a fraction, need to multiply it so that we don't get the same random number thousands of times.
    for nx = -math.floor(Nx/2), math.floor(Nx/2) do   --need to divied by two because we are including nx and -nx by using sines and cosines
       for ny = -math.floor(Ny/2), math.floor(Ny/2) do
         -- kx, ky = {2 pi nx/Lx, 2 pi ny/Ly}
         ksquared = math.pow(2*nx*math.pi /Lx,2) + math.pow(2*ny*math.pi/Ly,2)
         if ksquared > 0.0 then
            fIon = fIon + math.cos(2*math.pi*(nx*x/Lx + ny*y/Ly  + math.random() ))
          end  
       end
    end
    return {fIon}
 end

plasmaApp = Plasma.App {
    logToFile = true,
 
    tEnd        = 4000,          -- End time. RLW: I changed this, but didn't change nFrame below.  So we might get more frames.
    nFrame      = 400,             -- Number of output frames.  This is 0.5 frames in unit time --> 3 frames every Langmuir period. 
    nDistFuncFrame = 100,           -- Number of distribution function output frames 

    lower       = {0.0,0.0},             -- Configuration space lower left.
    upper       = {1.0,0.5},             -- Configuration space upper right.
    cells       = {32,16},               -- Configuration space cells.
    basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
    polyOrder   = pOrder,           -- Polynomial order.
    timeStepper = "rk3",            -- one of "rk2" or "rk3".
    cflFrac     = 0.9,
 
    -- Decomposition for configuration space.
    decompCuts = {32, 16},    -- Cuts in each configuration direction.
    useShared  = false,    -- If to use shared memory.
 
    -- Boundary conditions for configuration space.
    periodicDirs = {1,2}, -- Periodic directions.
    -- Integrated moment flag, compute quantities 1000 times in simulation.
    calcIntQuantEvery = 0.005,
    restartFrameEvery = 0.01, 
 
    -- Electrons.
    elc = Plasma.Species {
       charge = -1.0, mass = 1.0,
       -- Velocity space grid.
       lower = {-6.0*vtElc,-9.0*vtElc},
       upper = {20.0*vtElc,9.0*vtElc},
       cells = {416,144},
       decompCuts = {1,1},
       -- initial conditions
       init = function (t, xn)
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
          local fv = maxwellian2v({vx, vy}, vDriftElc, vtElc)
          --local dfv = sponEmissionSource({x,y,vx,vy},t, lx, nx, pOrder)[2]*maxwellian2v({vx, vy}, vDriftElc, vtElc)
          return fv
       end,
       evolve = true, -- Evolve species?
 
       diagnostics           = {"intM1i","intM2Thermal"},
       coll = Plasma.LBOCollisions{
        collideWith = {'elc'},
	    frequencies = {nuee},
      }
    },
 
    -- Ions.
    ion = Plasma.Species {
       charge = 1.0, mass = massRatio,
       -- Velocity space grid.
       lower = {-36.0*vtIon,-24.0*vtIon},
       upper = {72.0*vtIon,24.0*vtIon},
       cells = {54,24},
       decompCuts = {1,1},
       -- Initial conditions.
       init = function (t, xn)
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
          local fv = maxwellian2v({vx, vy}, vDriftIon, vtIon)
          local dfv = sponEmissionSource({x,y,vx,vy},t, lx, nx, pOrder)[1]*maxwellian2v({vx, vy}, vDriftIon, vtIon)
          return fv + perturbation*dfv    -- *(1.0 + perturbation*DrandKick({x,y},lx,nx,pOrder) ) RLW: I removed the noise
       end,
       source = Plasma.Source{
	 profile = function (t,xn)
         local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
         local dfv = sponEmissionSource({x,y,vx,vy},t, lx, nx, pOrder)[1]*maxwellian2v({vx,vy},{accion*t,0.0}, vtIon) 
         return noise*dfv
       end,        
       },

       evolve = true,    -- Evolve species?
 
       diagnostics           = {"intM1i","intM2Thermal"},
       coll = Plasma.LBOCollisions{
        collideWith = {'ion'},
	 frequencies = {nuii},
      }
    },

    -- Field solver.
    field = Plasma.Field {
       epsilon0 = 1.0,
       evolve = true, -- Evolve field?
       hasMagneticField = false,
    },

    externalField = Plasma.ExternalField {
      hasMagneticField = false,
      emFunc = function(t, xn)
         local extE_x, extE_y, extE_z = -0.00005, 0., 0.
         return extE_x, extE_y, extE_z
      end,
      evolve = false, -- Evolve field?
   }, 
 }
 -- Run application.
 plasmaApp:run()
