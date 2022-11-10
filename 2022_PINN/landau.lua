local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Logger = require "Logger"

local logger = Logger {logToFile = True}
local log = function(...)
   logger(string.format(...))
   logger("\n")
end

-- parameters
-- background values before perturbation
local q_e, q_i = -1, 1
local m_e, m_i = 1, 1836
local n_e, n_i = 1, 1
local T_e, T_i = 1, 1 -- temperature
local vth_e, vth_i = math.sqrt(T_e / m_e), math.sqrt(T_i / m_i) -- thermal speed

local c = 10 * vth_e -- artificially-decreased speed of light
local epsilon_0 = 1
local mu_0 = 1 / (epsilon_0 * c ^ 2)

-- derived parameters
local omega_pe = math.sqrt((n_e * q_e ^ 2) / (epsilon_0 * m_e)) -- plasma freq.
local lambda_D = math.sqrt((epsilon_0 * T_e) / (n_e * q_e ^ 2)) -- Debye length
local Lx = 2 * math.pi / 0.6 * lambda_D -- simulation domain length

-- perturbations
local lambda1 = Lx / 1 -- wavelength of mode 1
local lambda2 = Lx / 2 -- wavelength of mode 2
local kx1 = 2 * math.pi / lambda1 -- wavenumber of mode 1
local kx2 = 2 * math.pi / lambda2 -- wavenumber of mode 2
-- it important to make sure the density fluctuation amplitudes are much
-- smaller than 1 to mak
local A1 = 5e-2 -- relative density fluctuation magnitude of mode 1
local A2 = 4 * A1 / kx1 * kx2 -- relative density fluctuation magnitude of mode 2
local phi2 = 0.38716 -- a random phase shift of mode 2

local tEnd = 30 / omega_pe
-- set maximumDt slightly smaller than dt given by the cfl condition so we can
-- safely use a constant dt during the simulation
local maximumDt = 0.001
local nFrame = tEnd / maximumDt + 1 -- output every step
nFrame = 10

log("%30s = %g", "Lx", Lx)
log("fluctuation wavenumber and wavelengths")
log("%30s = %g", "kx1", kx1)
log("%30s = %g", "kx2", kx2)
log("%30s = %g", "lambda1", lambda1)
log("%30s = %g", "lambda2", lambda2)
log("density field fluctuation relative magnitudes:")
log("%30s = %g", "A1", A1)
log("%30s = %g", "A2", A2)
log("E field fluctuation relative magnitudes:")
log("%30s = %g", "A1/kx1", A1 / kx1)
log("%30s = %g", "A2/kx2", A2 / kx2)

-- Maxwellian distribution function
local function maxwellian1D(n, vth, vx)
   return n / math.sqrt(2 * math.pi * vth * vth) *
              math.exp(-(vx) ^ 2 / (2 * vth * vth))
end

-- create the simulation application
local sim = Plasma.App {
   logToFile = false,

   tEnd = tEnd, -- simulation ending time
   nFrame = nFrame, -- number of time frames to ouput
   maximumDt = maximumDt,
   lower = {0}, -- configuration space lower corner in each direction
   upper = {Lx}, -- configuration space upper corner in each direction
   cells = {128}, -- configuration space number of cells in each direction
   decompCuts = {4}, -- number of processors along each direction
   basis = "serendipity", -- polynomial basis type for the DG algorithm
   polyOrder = 2, -- polynomial basis order for the DG algorithm
   timeStepper = "rk3", -- time integrator type; one of "rk2" or "rk3"
   cflFrac = 1, -- CFL " fraction "; usually 1
   periodicDirs = {1}, -- configuration space directions w. periodic boundaries

   -- eclectrons
   elc = Plasma.Species {
      charge = q_e,
      mass = m_e,
      -- velocity space grid
      lower = {-6 * vth_e},
      upper = {6 * vth_e},
      cells = {32},
      -- initial conditions
      init = function(t, xn)
         local x, vx = xn[1], xn[2]
         local n = n_e
         n = n + n_e * A1 * math.cos(kx1 * x)
         n = n + n_e * A2 * math.cos(kx2 * x + phi2)
         return maxwellian1D(n, vth_e, vx)
      end,
      evolve = true,
      diagnostics = {"M0", "M1i", "M2ij", "M2", "M3i"}
   },

   -- ions (immobile background)
   ion = Plasma.Species {
      charge = q_i,
      mass = m_i,
      -- velocity space grid
      lower = {-6 * vth_i},
      upper = {6 * vth_i},
      cells = {32},
      -- initial conditions
      init = function(t, xn)
         local x, vx = xn[1], xn[2]
         return maxwellian1D(n_i, vth_i, vx)
      end,
      evolve = false,
      diagnostics = {"M0", "M1i", "M2ij", "M2", "M3i"}
   },

   -- electromagnetic field
   field = Plasma.Field {
      epsilon0 = epsilon_0,
      mu0 = mu_0,
      hasMagneticField = false,
      evolve = true,
   }
}

-- run the simulation
sim:run()

