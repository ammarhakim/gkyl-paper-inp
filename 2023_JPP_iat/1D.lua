-- Gkyl ------------------------------------------------------------------------
-- Z.Liu 6/25/2022
-- Mass=25, Temp=50
-- E1, E_ext = 1.0e-5
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

-- electron parameters
vDriftElc = 0.0   ---modified from 0.159
vtElc = 0.02
-- ion parameters
vDriftIon = 0.0
vtIon = 0.0002828      ---- modified from 0.001 (use Te = 50 Ti)
tempRatio = 50.0
-- mass ratio
massRatio = 100.0  ----modified from 25 

perturbation = 1.0e-6 -- distribution function perturbation
noise = 1.0e-6

omegaPe = 1.0
omegaPi = omegaPe/math.sqrt(massRatio)
lambdaDe = vtElc / omegaPe

-- Collision frequencies.
nuee = 0.0001
nuei = 0.0 --nuee              --RLW 
nuii = 2.0e-5  --should have been nuee/math.sqrt(massRatio)*math.pow(tempRatio,1.5)
nuie = 0.0

lx = 1.0
nx = 64

pOrder = 2

accion = -5.0e-7

local function maxwellian1v(v, vDrift, vt)
    if v-vDrift > 4.0*vtElc then	
       return 0.00001*1.0/math.sqrt(2*math.pi*vt^2)*math.exp(-((v-vDrift)^2)/(2*vt^2))
    elseif v-vDrift <- 4.0*vtElc then
       return 0.00001*1.0/math.sqrt(2*math.pi*vt^2)*math.exp(-((v-vDrift)^2)/(2*vt^2))
    else
       return 1.0/math.sqrt(2*math.pi*vt^2)*math.exp(-((v-vDrift)^2)/(2*vt^2))
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
    local Nx = ncells_table*(p+1) -- Number of spatial degrees of freedom along x
    local Lx = lx_table
    local x, vx = x_table[1], x_table[2]
    local fIon = 0.0
    local ksquared = 0.0
    math.randomseed(math.floor(1000000*t)) --I want all xs and ys (and vx, vy) to see the same random phase at a given time step.  Since t will be a fraction, need to multiply it so that we don't get the same random number thousands of times.
    for nx = -math.floor(Nx/2), math.floor(Nx/2) do   --need to divied by two because we are including nx and -nx by using sines and cosines
         fIon = fIon + math.cos(2*math.pi*(nx*x/Lx  + math.random() ))
      end
    return {fIon}
 end

plasmaApp = Plasma.App {
    logToFile = true,
 
    tEnd        = 10000,          -- End time. RLW: I changed this, but didn't change nFrame below.  So we might get more frames.
    nFrame      = 1000,             -- Number of output frames.  This is 0.5 frames in unit time --> 3 frames every Langmuir period. 
    nDistFuncFrame = 100,           -- Number of distribution function output frames 

    lower       = {0.0},             -- Configuration space lower left.
    upper       = {lx},             -- Configuration space upper right.
    cells       = {nx},               -- Configuration space cells.
    basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
    polyOrder   = pOrder,           -- Polynomial order.
    timeStepper = "rk3",            -- one of "rk2" or "rk3".
    cflFrac     = 0.9,
 
    -- Decomposition for configuration space.
    decompCuts = {32},    -- Cuts in each configuration direction.
    useShared  = false,    -- If to use shared memory.
 
    -- Boundary conditions for configuration space.
    periodicDirs = {1}, -- Periodic directions.
    -- Integrated moment flag, compute quantities 1000 times in simulation.
 
    -- Electrons.
    elc = Plasma.Species {
       charge = -1.0, mass = 1.0,
       -- Velocity space grid.
       lower = {-6.0*vtElc},
       upper = {18.0*vtElc},
       cells = {768},
       decompCuts = {1},
       -- initial conditions
       init = function (t, xn)
          local x, vx = xn[1], xn[2]
          local fv = maxwellian1v(vx, vDriftElc, vtElc)
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
       lower = {-48.0*vtIon},
       upper = {72.0*vtIon},
       cells = {240},
       decompCuts = {1},
       -- Initial conditions.
       init = function (t, xn)
          local x, vx = xn[1], xn[2]
          local fv = maxwellian1v(vx, vDriftIon, vtIon)
          local dfv = sponEmissionSource({x,vx},t, lx, nx, pOrder)[1]*maxwellian1v(vx, vDriftIon, vtIon)
          return fv + perturbation*dfv    -- *(1.0 + perturbation*DrandKick({x,y},lx,nx,pOrder) ) RLW: I removed the noise
       end,
       source = Plasma.Source{
	 profile = function (t,xn)
         local x, vx = xn[1], xn[2]
         local dfv = sponEmissionSource({x,vx},t, lx, nx, pOrder)[1]*maxwellian1v(vx,accion*t, vtIon) 
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
