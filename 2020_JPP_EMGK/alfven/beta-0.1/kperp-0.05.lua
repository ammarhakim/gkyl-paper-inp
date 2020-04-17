-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

TiTe = 1.0
knumber = 1.0
beta = 0.1
kperp2 = 0.05^2
alpha = 1e-6
Ti0 = 1.0
Te0 = Ti0/TiTe
n0 = 1.0
ni0, ne0 = n0, n0
B0 = 1
omega = 1.8
--print("expected freq = ", omega)
--print("expected period = ", 2*math.pi/omega)

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 2*math.pi/omega*5, --5/knumber, -- end time
   nFrame = 1, -- number of output frames
   lower = {-math.pi/knumber}, -- configuration space lower left
   upper = {math.pi/knumber}, -- configuration space upper right
   cells = {32}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {32}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = -1.0,
      mass = 1.0,
      -- velocity space grid
      lower = {-8.0},
      upper = {8.0},
      cells = {64},
      decompCuts = {1},
      -- initial conditions
      initBackground = {"maxwellian",
              density = function (t, xn)
                 return ne0
              end,
              temperature = function (t, xn)
                 return Te0
              end,
             },
      init = {"maxwellian",
              density = function (t, xn)
	         local x, v = xn[1], xn[2]
                 local k = knumber
                 return ne0*(1 + alpha*math.cos(k*x))
              end,
              temperature = function (t, xn)
                 return Te0
              end,
             },
      evolve = true, -- evolve species?
   },

   stationaryIon = Plasma.GkSpecies {
      charge = 1.0,
      mass = 1.0,
      -- velocity space grid
      lower = {-8.0},
      upper = {8.0},
      cells = {64},
      decompCuts = {1},
      init = {"maxwellian",
              density = function (t, xn)
	         local x, v = xn[1], xn[2]
                 local k = knumber
                 return ni0
              end,
              temperature = function (t, xn)
                 return Ti0
              end,
             },
      evolve = false, -- evolve species?
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve field?
      isElectromagnetic = true,
      kperp2 = kperp2,
      mu0 = beta,
      polarizationWeight = 1.0,
      discontinuousPhi = false,
      discontinuousApar = true,
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         return B0
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
