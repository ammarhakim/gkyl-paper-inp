-- Gkyl ------------------------------------------------------------------------
-- 3D test simulation of LAPD
--------------------------------------------------
-- Load magnetic field data
local lines = io.lines("High_vA_profile.txt") -- open file as lines
local LAPDTable = {} -- new table with columns and rows as tables[n_column][n_row]=value
for line in lines do -- row iterator
    local i = 1 -- first column
    for value in (string.gmatch(line, "[^%s]+")) do  -- tab separated values
--  for value in (string.gmatch(line, '%d[%d.]*')) do -- comma separated values
        LAPDTable[i]=LAPDTable[i]or{} -- if not column then create new one
        LAPDTable[i][#LAPDTable[i]+1]=tonumber(value) -- adding row value
        i=i+1 -- column iterator
    end
end

local zLAPD = LAPDTable[1]
local BzLAPD = LAPDTable[2] -- In units of kGauss! Radial coords 2-11
local offsetBz = 1
-- End load data
-- Functions for finding nearest indicies in data and interpolating
local findNearestIndex = function(x, x0)
      local nx = #x
      local idxLo = 1
      while (x[idxLo+1] <= x0 and idxLo+1 < nx) do
      	    idxLo = idxLo + 1
      end
      local idxUp = idxLo + 1
      return idxLo, idxUp
end
local linearInterp = function(f, x, x0)
      local idxLo, idxUp = findNearestIndex(x, x0)
      local x0l = x0
      if (x0 > x[#x]) then
         x0l = x[#x]
      end
      local y = f[idxLo] + (x0l - x[idxLo]) * (f[idxUp] - f[idxLo]) / (x[idxUp] - x[idxLo])
      return y
end
local biLinearInterp = function(f, x, x0, y, y0, offset)
      local idxLoX, idxUpX = findNearestIndex(x, x0)
      local idxLoY, idxUpY = findNearestIndex(y, y0)
      idxLoYTab = idxLoY + offset -- deal with 2D LAPDTable indexing
      idxUpYTab = idxUpY + offset
      local x0l, y0l = x0, y0
      if (x0 > x[#x]) then
      	 x0l = x[#x]
      end
      if (y0 > y[#y]) then
      	 y0l = y[#y]
      end
      local f0 = f[idxLoYTab][idxLoX] + (x0l - x[idxLoX]) * (f[idxLoYTab][idxUpX] - f[idxLoYTab][idxLoX]) / (x[idxUpX] - x[idxLoX])
      local f1 = f[idxUpYTab][idxLoX] + (x0l - x[idxLoX]) * (f[idxUpYTab][idxUpX] - f[idxUpYTab][idxLoX]) / (x[idxUpX] - x[idxLoX])

      local y = f0 + (y0l - y[idxLoY]) *  (f1 - f0) / (y[idxUpY] - y[idxLoY])
      return y
end
-- example Bz = biLinearInterp(LAPDTable, zLAPD, 0., rLAPD, 0.2, offsetBz)
----------------------------------------------

-- local Moments = require "Moments"
local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local Constants = require "Lib.Constants"
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = true
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end


-- physical parameters
gasGamma = 5./3.
local cFac = 1000
-- Universal constant parameters.
local eps0, eV = Constants.EPSILON0*cFac, Constants.ELEMENTARY_CHARGE
local mu0      = Constants.MU0
local m_e, m_p  = Constants.ELECTRON_MASS, Constants.PROTON_MASS

lightSpeed = 1/math.sqrt(mu0*eps0) -- speed of light

Te_Ti = 5.0 -- ratio of electron to ion temperature
elcTemp = 7.*eV
ionTemp = elcTemp / Te_Ti
n0 = 7.0e18 -- initial number density

ionMass = 4*m_p -- ion mass, Helium
ionCharge = 1*eV -- ion charge, singly ionized He
elcMass = ionMass / 100. -- m_e -- electron mass
elcCharge = -eV -- electron charge

B0 =  linearInterp(BzLAPD, zLAPD, 0.) / 10000. -- convert to Tesla
P0 = B0*B0 / (2*mu0) + n0*(elcTemp + ionTemp) -- pressure at antenna end

-- Alfven speeds
vAIon = B0/math.sqrt(mu0*n0*ionMass)
vAElc = B0/math.sqrt(mu0*n0*elcMass)

-- Thermal speeds
vtElc = math.sqrt(2.*elcTemp / elcMass) --electron thermal velocity
vtIon = math.sqrt(2.*ionTemp / ionMass) --ion thermal velocity

-- Ion beta
beta = vtIon^2 / vAIon^2

-- Elc beta
betaElc = vtElc^2 / vAElc^2

-- Beta bar
betaBar = vtElc^2 / vAIon^2

-- plasma and cyclotron frequencies
wpe = math.sqrt(ionCharge^2*n0/(eps0*elcMass))
wpi = math.sqrt(ionCharge^2*n0/(eps0*ionMass))
omegaCe = ionCharge*B0/elcMass
omegaCi = ionCharge*B0/ionMass

-- inertial length
de = vAElc/omegaCe
di = vAIon/omegaCi

-- gyroradii
rhoi = vtIon/omegaCi
rhoe = vtElc/omegaCe
rhos = math.sqrt(elcTemp/ionMass) / omegaCi

-- antenna params
J0 = 1.0e4   -- Amps/m^3.
driveFreq = 7.64e4
antRamp = 0.25/driveFreq
tAntOff = 1.5/driveFreq - antRamp
lAnt = 0.2178
kAntx = 2*math.pi/lAnt


-- domain size and simulation time
Lr = 0.6
LzSt = -10.
LzEnd = 18.
Lz = LzEnd - LzSt

nr = 64
nz = 700

ncx = 8
ncy = 8
ncz = 24

tEnd = 150.0/omegaCi
nFrames = 150

dz = Lz / nz
zFirstEdge = dz

tTransit = LzEnd / vAIon

log("%50s = %g", "n0", n0)
log("%50s = %g", "mi/me", ionMass / elcMass)
log("%50s = %g", "wpe/OmegaCe", wpe / omegaCe)
log("%50s = %g", "proton beta", beta)
log("%50s = %g", "electron beta", betaElc)
log("%50s = %g", "beta_bar", betaBar)
log("%50s = %g", "vAIon", vAIon)
log("%50s = %g", "Alfven transit time in OmegaCi", tTransit*omegaCi)
log("%50s = %g", "c", lightSpeed)
log("%50s = %g", "vte/c", vtElc / lightSpeed)
log("%50s = %g", "vti/c", vtIon / lightSpeed)
log("%50s = %g", "electron plasma frequency (wpe) ", wpe)
log("%50s = %g", "electron cyclotron frequency (OmegaCe) ", omegaCe)
log("%50s = %g", "ion plasma frequency (wpi) ", wpi)
log("%50s = %g", "ion cyclotron frequency (OmegaCi) ", omegaCi)
log("%50s = %g", "electron inertial length (de) ", de)
log("%50s = %g", "ion inertial length (di) ", di)
log("%50s = %g", "electron gyroradius (rhoe) ", rhoe)
log("%50s = %g", "ion gyroradius (rhoi) ", rhoi)
log("%50s = %g", "ion sound gyroradius (rhos) ", rhos)
log("%50s = %g", "antenna freqeucny (Hz) ", driveFreq)
log("%50s = %g", "antenna perpendicular wavelength (m) ", 2*math.pi/kAntx)
log("%50s = %g", "Number of grid cells per di in z", nz/(Lz/di))
log("%50s = %g", "Number of grid cells per de in z", nz/(Lz/de))
log("%50s = %g", "Number of grid cells per di in x", nr/(2*Lr/di))
log("%50s = %g", "Number of grid cells per de in x", nr/(2*Lr/de))
log("%50s = %g", "tFinal ", tEnd)
log("%50s = %g", "End time in inverse ion cyclotron periods", tEnd*omegaCi)

momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrames,
   lower = {-Lr, -Lr, LzSt}, -- configuration space lower left
   upper = {Lr, Lr, LzEnd}, -- configuration space upper right
   cells = {nr, nr, nz}, -- configuration space cells
   timeStepper = "fvDimSplit", 
   -- decomposition for configuration space
   decompCuts = {ncx, ncy, ncz}, -- cuts in each configuration direction

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions

   -- electrons
   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,

      equation    = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      -- initial conditions
      init = function (t, xn)
	 local x, y, z = xn[1], xn[2], xn[3]
	 local r = math.sqrt(x*x + y*y)
	 local A = 1.
	 local ne = n0*A
	 local Te = elcTemp*A

	 local vdrift_x = 0.
         local vdrift_y = 0.
         local vdrift_z = 0.

	 local rhoe = ne*elcMass
         local exmom = vdrift_x*rhoe
	 local eymom = vdrift_y*rhoe
         local ezmom = vdrift_z*rhoe
	 local ere = ne*Te/(gasGamma-1) + 0.5*exmom*exmom/rhoe + 0.5*eymom*eymom/rhoe + 0.5*ezmom*ezmom/rhoe
	 
	 return rhoe, exmom, eymom, ezmom, ere
      end,
      evolve = true, -- Evolve species?
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcy = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcz = { Moments.Species.bcCopy, Moments.Species.bcCopy },
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
	 local r = math.sqrt(x*x + y*y)
	 local A = 1.
	 local ni = n0*A
	 local Ti = ionTemp*A

       	 local vdrift_x = 0.
         local vdrift_y = 0.
	 local vdrift_z = 0.

	 local rhoi = ni*ionMass
         local ixmom = vdrift_x*rhoi
	 local iymom = vdrift_y*rhoi
         local izmom = vdrift_z*rhoi
	 local eri = ni*Ti/(gasGamma-1) + 0.5*ixmom*ixmom/rhoi + 0.5*iymom*iymom/rhoi + 0.5*izmom*izmom/rhoi
	 
	 return rhoi, ixmom, iymom, izmom, eri

      end,
      evolve = true, -- Evolve species?
      bcx = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcy = { Moments.Species.bcWall, Moments.Species.bcWall },
      bcz = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   field = Moments.Field {
      epsilon0 = eps0, mu0 = mu0,
      init = function (t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
	 
         return  0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,

      evolve = true, -- Evolve species?
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },
      bcy = { Moments.Field.bcReflect, Moments.Field.bcReflect },
      bcz = { Moments.Field.bcReflect, Moments.Field.bcReflect },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "time-centered",
      hasStaticField = true,
      staticEmFunction = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         
         local Bz = linearInterp(BzLAPD, zLAPD, z) / 10000.
	 if z < 0. then
		Bz = B0
	 end
 	 return 0.0, 0.0, 0.0, 0.0, 0.0, Bz
      end,

       -- Additional source terms.
      -- Not enabled here; for demo purpose.
      -- Note: Dp are c arrays and use 0-based indices; xc and qbym are lua
      --       arrays and use 1-based indices
      hasAuxSourceFunction = true,
      auxSourceFunction    = function (self, xc, t, epsilon0, qbym, fDp, emDp, auxSrcDp)
         local x, y, z = xc[1], xc[2], xc[3]
         local nFluids = #qbym
	 local r = math.sqrt(x*x + y*y)
        
         -- Auxiliary source for currents.
         for s=0, 0 do
            auxSrcDp[s*3+0] = 0
            auxSrcDp[s*3+1] = 0
            auxSrcDp[s*3+2] = 0
         end

         -- Auxiliary source for E field.
         auxSrcDp[nFluids*3+0] = 0
         auxSrcDp[nFluids*3+1] = 0
         auxSrcDp[nFluids*3+2] = 0
	 
         
         if (z < dz/2. and z >= -dz/2.) then
		  antTemp = -J0*math.sin(2.*math.pi*driveFreq*t)*math.sin(.5*math.pi*math.min(1,t/antRamp))^2.*math.cos(0.5*math.pi*math.min(1, math.max(0,(t-tAntOff)/antRamp)))^2.
        	  antSpat =  math.sin(kAntx*x)
		  antSpatShape = 0.5*(1 - math.tanh((math.abs(y) - lAnt/2.)*10/lAnt))*0.5*(1 - math.tanh((math.abs(x) - 2*lAnt/3.)*10/lAnt))
		  auxSrcDp[nFluids*3+2] = antTemp*antSpat*antSpatShape/eps0
         end
      end,


   },   
}
-- run application
momentApp:run()
