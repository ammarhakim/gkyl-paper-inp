-- Gkyl ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"

-- Constants
elcCharge = -1.0
ionCharge = 1.0
elcMass = 1.0
ionMass = 36.0
Te_Ti = 0.1

-- Initial conditions
n0 = 1.0
nbOverN0 = 0.001
vtElc = 0.06
elcTemp = vtElc*vtElc*elcMass/2.0
ionTemp = elcTemp/Te_Ti
vtIon = math.sqrt(2.0*ionTemp/ionMass)

-- electron beta
beta = 1.0/11.0 --total beta = 1.0 = ion beta + electron beta
vAe = vtElc/math.sqrt(beta)
B0 = vAe --derived from normalization of elcMass and n0
vAi = vAe/math.sqrt(ionMass)

omegaCi = ionCharge*B0/ionMass
omegaCe = ionCharge*B0/elcMass

larmor_elc = vtElc/omegaCe
larmor_ion = vtIon/omegaCi

-- perturbation
noiseAmp = 0.0001
mode = 8.0

l = larmor_ion -- Current sheet width
Lx = 6.4*l
Ly = 12.8*l
Nx = 128
Ny = 256
endTime = 8.0/omegaCi

-- Maxwellian in 2x2v
local function maxwellian2D(n, vx, vy, ux, uy, temp, mass)
   local v2 = (vx - ux)^2 + (vy - uy)^2
   return n/(2.0*math.pi*temp/mass)*math.exp(-mass*v2/(2.0*temp))
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = endTime, -- end time
   nFrame = 40, -- number of output frames
   lower = {-Lx/2, -Ly/2}, -- configuration space lower left
   upper = {Lx/2, Ly/2}, -- configuration space upper right
   cells = {Nx, Ny}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {32,32}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1}, -- periodic directions

   -- electrons 
   elc = Plasma.VlasovSpecies {
      charge = elcCharge, mass = elcMass,
      -- velocity space grid
      lower = {-8.0*vtElc, -8.0*vtElc},
      upper = {8.0*vtElc, 8.0*vtElc},
      cells = {32, 32},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
       local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
       local cosh = math.cosh
       local sin = math.sin
       local cos = math.cos
       local sech2 = (1.0/cosh(y/l))^2
       local Pi = math.pi
       local _2pi = 2.0*Pi

       local n = n0*sech2
       local nBackground = n0*nbOverN0
       local TeFrac = elcTemp / (elcTemp + ionTemp)
       local TiFrac = 1.0 - TeFrac
       local i_y = 1.0
       local i_x = mode
       local JxNoise = -(i_y*Pi/Ly)*sin(i_y*Pi*y/Ly)*sin(i_x*2.0*Pi*x/Lx)
       local JyNoise = -(i_x*2.0*Pi/Lx)*cos(i_y*Pi*y/Ly)*cos(i_x*2.0*Pi*x/Lx)

       JxNoise = noiseAmp*JxNoise/mode
       JyNoise = noiseAmp*JyNoise/mode

       local Jx  = (B0/l)*(-sech2)  + JxNoise
       local Jy  = JyNoise
       local uxe = Jx*TeFrac/(elcCharge*n)
       local uxi = Jx*TiFrac/(ionCharge*n)
       local uye = Jy*TeFrac/(elcCharge*n)
       local uyi = Jy*TiFrac/(ionCharge*n)
       local fv= maxwellian2D(n, vx, vy, uxe, uye, elcTemp, elcMass)+maxwellian2D(nBackground, vx, vy, 0.0, 0.0, elcTemp, elcMass)
       return fv
      end,
      -- boundary conditions
      bcy = { Plasma.VlasovSpecies.bcReflect, Plasma.VlasovSpecies.bcReflect },
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },
   -- protons
   ion = Plasma.VlasovSpecies {
      charge = ionCharge, mass = ionMass,
      -- velocity space grid
      lower = {-6.0*vtIon, -6.0*vtIon},
      upper = {6.0*vtIon, 6.0*vtIon},
      cells = {24, 24},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn)
       local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
       local cosh = math.cosh
       local sin = math.sin
       local cos = math.cos
       local sech2 = (1.0/cosh(y/l))^2
       local Pi = math.pi
       local _2pi = 2.0*Pi

       local n = n0*sech2
       local nBackground = n0*nbOverN0
       local TeFrac = elcTemp / (elcTemp + ionTemp)
       local TiFrac = 1.0 - TeFrac
       local i_y = 1.0
       local i_x = mode
       local JxNoise = -(i_y*Pi/Ly)*sin(i_y*Pi*y/Ly)*sin(i_x*2.0*Pi*x/Lx)
       local JyNoise = -(i_x*2.0*Pi/Lx)*cos(i_y*Pi*y/Ly)*cos(i_x*2.0*Pi*x/Lx)

       JxNoise = noiseAmp*JxNoise/mode
       JyNoise = noiseAmp*JyNoise/mode

       local Jx  = (B0/l)*(-sech2)  + JxNoise
       local Jy  = JyNoise
       local uxe = Jx*TeFrac/(elcCharge*n)
       local uxi = Jx*TiFrac/(ionCharge*n)
       local uye = Jy*TeFrac/(elcCharge*n)
       local uyi = Jy*TiFrac/(ionCharge*n)
       local fv= maxwellian2D(n, vx, vy, uxi, uyi, ionTemp, ionMass)+maxwellian2D(nBackground, vx, vy, 0.0, 0.0, ionTemp, ionMass)
       return fv
      end,
      -- boundary conditions
      bcy = { Plasma.VlasovSpecies.bcReflect, Plasma.VlasovSpecies.bcReflect },
      evolve = true, -- evolve species?
      diagnosticMoments = { "M0", "M1i", "M2", "M2ij", "M3i" },
      diagnosticIntegratedMoments = { "intM0", "intM1i", "intM2Flow", "intM2Thermal" },
   },

   -- field solver
   field = Plasma.MaxwellField {
      epsilon0 = 1.0, mu0 = 1.0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local tanh = math.tanh
         local sin = math.sin
         local cos = math.cos
         local Pi = math.pi
         local _2pi = 2.0*Pi

         local i_y = 1.0
         local i_x = mode
         local BzNoise = cos(i_y*Pi*y/Ly)*sin(i_x*2.0*Pi*x/Lx)
         BzNoise = noiseAmp*BzNoise/mode
         local Bzb = -B0*tanh(y/l)
         local Bz = Bzb + BzNoise
	 return 0.0, 0.0, 0.0, 0.0, 0.0, Bz
      end,
      bcy = { Plasma.MaxwellField.bcReflect, Plasma.MaxwellField.bcReflect },
      evolve = true, -- evolve field?
   },
}
-- run application
plasmaApp:run()
