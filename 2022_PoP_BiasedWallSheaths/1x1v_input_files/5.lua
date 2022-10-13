-- Gkyl --------------------------------------------------------------
-- Basic sheath simulation -------------------------------------------
local Plasma = (require "App.PlasmaOnCartGrid").VlasovMaxwell()
local const = require "Lib.Constants"

-- Parameters we want to change in these input files
local V_left = 5.e3     -- in V   (though you will typically want to use ~kV values)
local Robertson_profile_name = "5"

-- Set constants
local epsilon_0 = const.EPSILON0
local mu_0 = const.MU0
local q_e, q_i = -const.ELEMENTARY_CHARGE, const.ELEMENTARY_CHARGE
local m_e, m_i = const.ELECTRON_MASS, const.MASS_UNIT  -- Using proton-electron plasma

-- Set densities and temperatures based on FuZe parameters
-- See Zhang et al (2019) PRL
-- https://doi.org/10.1103/PhysRevLett.122.135001
local n_e, n_i = 1.1e23, 1.1e23   --m^-3
local T_e, T_i = 2000.0*const.ELEMENTARY_CHARGE, 2000.0*const.ELEMENTARY_CHARGE    -- in J, 2 keV

-- Calculate thermal speeds, plasma frequency, and Debye length
local vth_e, vth_i = math.sqrt(T_e/m_e), math.sqrt(T_i/m_i)
local omega_pe = math.sqrt((n_e * q_e^2)/(epsilon_0*m_e))
local lambda_D = math.sqrt((epsilon_0 * T_e)/(n_e * q_e^2))

-- Calculate Bohm speed 
-- This is used for Robertson profile initial conditions (more on this lower down)
-- Note that we don't include T_i in this calculation because Robertson assumes T_e >> T_i
local cs = math.sqrt((T_e)/m_i)

-- Initialize drift velocities to zero
local vd_e, vd_i = 0.0 , 0.0

-- Set grid parameters
local domainLength_norm = 256
local nx = 2
local domainLength = domainLength_norm*lambda_D

-- Set the "length" of the sheath region
-- This is the region where we will NOT reintroduce plasma via the source function
local sheathLength = 28.0*lambda_D

-- Set the length of the linearly increasing and decreasing portions of the source region
local linearLength = 100.0*lambda_D

-- Set the length of the center constant portion of the source region (this should be be equal to domainLength-2*sheathLength-2*linearLength)
local centerLength = domainLength - 2*(sheathLength + linearLength)


-- Make a function for making Maxwellians
local function maxwellian(n, u, vth, v)
   return n/math.sqrt(2*math.pi*vth*vth)*math.exp(-(v-u)^2/(2*vth*vth))
end

-- Make functions to calculate the linearly increasing and decreasong part of source term (used to maintain particle conservation in domain)
-- The integral of the entire source function integrated over x should yield 1
-- Note that there will be an additional long region of constant source (profile is trapezoidal)
-- To make these function, first, define a slope
local sourceSlope = 1.0/(linearLength*(linearLength+centerLength))
local function fL(x)
   return sourceSlope*(x-sheathLength)
end

local function fR(x)
   return -sourceSlope*(x-sheathLength-linearLength-centerLength)+1.0/(linearLength+centerLength)
end

-- Set collision frequency values
local mfp = 50.*lambda_D
local nu_ee = vth_e/mfp
local nu_ei = nu_ee
local nu_ii = vth_i/mfp
local nu_ie = (m_e/m_i)*nu_ee

-- Make function to set up the collision profiles
-- Collisional profile is needed to have a collisional pre-sheath and a collisionless sheath
-- This profile is then multiplied by a scalar collision frequency

-- First make the base function that is modeled after Kolter's function
local function fk(x)
   return 1.0/(1.0 + math.exp(x/(12.0*lambda_D) - 8.0/1.5))
end

local function collProfile(x)
   return fk(-x+128.0*lambda_D) + fk(x-domainLength+128.0*lambda_D) - 1
end

-- Make collision profile functions for each type of collision we are considering
local function nu_eeProfile(t, xn)
   local x = xn[1]
   return nu_ee*collProfile(x)
end
local function nu_eiProfile(t, xn)
   local x = xn[1]
   return nu_ei*collProfile(x)
end
local function nu_iiProfile(t, xn)
   local x = xn[1]
   return nu_ii*collProfile(x)
end
local function nu_ieProfile(t, xn)
   local x = xn[1]
   return nu_ie*collProfile(x)
end

-- Load data for initial left and right Robertson profiles
-- All Robertson profiles calculated using method laid out in Robertson (2013) PPCF
-- dx.doi.org/10.1088/0741-3335/55/9/093001
-- Note that normalized ion particle flux at the wall is set to 0.55 (based on initial simulation results)
-- The wall potential for Robertson method are calculated using the biased wall potential formulas found in 
-- Stangeby (2000) textbook, Plasma Boundary of Magnetic Sheaths, Section 2.6
-- These profiles are used as initial conditions as they provide better approximations than a uniform plasma (i.e. reach steady state more quickly)
local fh = io.open("xL_" .. Robertson_profile_name .. ".txt")
local xL, i = {}, 1
for l in fh:lines() do
   xL[i] = l*1.0
   i = i+1
end

local fh = io.open("niL_" .. Robertson_profile_name .. ".txt")
local n_i0L, i = {}, 1
for l in fh:lines() do
   n_i0L[i] = l*1.0
   i = i+1
end

local fh = io.open("neL_" .. Robertson_profile_name .. ".txt")
local n_e0L, i = {}, 1
for l in fh:lines() do
   n_e0L[i] = l*1.0
   i = i+1
end

local fh = io.open("uiL_" .. Robertson_profile_name .. ".txt")
local uiL, i = {}, 1
for l in fh:lines() do
   uiL[i] = l*1.0
   i = i+1
end

-- Get the total number of elements for the left profile arrays
local numL = i-1

local fh = io.open("xR_" .. Robertson_profile_name .. ".txt")
local xR, i = {}, 1
for l in fh:lines() do
   xR[i] = l*1.0
   i = i+1
end

local fh = io.open("niR_" .. Robertson_profile_name .. ".txt")
local n_i0R, i = {}, 1
for l in fh:lines() do
   n_i0R[i] = l*1.0
   i = i+1
end

local fh = io.open("neR_" .. Robertson_profile_name .. ".txt")
local n_e0R, i = {}, 1
for l in fh:lines() do
   n_e0R[i] = l*1.0
   i = i+1
end

local fh = io.open("uiR_" .. Robertson_profile_name .. ".txt")
local uiR, i = {}, 1
for l in fh:lines() do
   uiR[i] = l*1.0
   i = i+1
end

-- Get the total number of elements for the right profile arrays
local numR = i-1

-- Get the maximum lengths for the left and right profiles
local xL_max = xL[numL]
local xR_max = xR[numR]

sim = Plasma.App {
   logToFile = false,

   tEnd = 20000./omega_pe, -- end time
   nFrame = 20, -- number of output frames
   lower = {0.0*lambda_D}, -- configuration space lower left
   upper = {domainLength}, -- configuration space upper right
   cells = {nx*domainLength_norm}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 2, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"

   -- decomposition for configuration space
   decompCuts = {128}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {}, -- periodic directions
   calcIntQuantEvery = 0.05,
   -- electrons
   elc = Plasma.Species {
      charge = q_e, mass = m_e,
      -- velocity space grid
      lower = {-6.0*vth_e},
      upper = {6.0*vth_e},
      cells = {64},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                -- get position normalized by Debye length 
                local x = xn[1]/lambda_D      
                
                -- Make an index for where we want to take the values from
                local val_index = 1
                
                -- Split this into three different regions, left, center, and right
                if x < xL_max then
                   -- We are in the left region                
                   -- Iterate through xL to determine where we should get the data from
                   for i = 1,numL do 
                      -- We want to use the value that corresponds to the highest value of xL below x
                      -- One error that comes up with this method is that some values of x initially might be negative
                      -- To fix this, we will set up an error catching early on to set and negative x value points to the first index
                      -- This is accomplished by skipping any x <= 0                                              
                      if xL[i] > x and x > 0 then
                         val_index = i-1                         
                         break
                      end
                   end
                   return n_e0L[val_index]*n_e  
                elseif x > domainLength_norm-xR_max+xR[1] then
                  -- We are in the right region
                  -- Note that we had to add by xR[1] since the first value is not at 0 for the right profiles
                  for i = 1,numR do 
                      -- We want to use the value that corresponds to the highest value of xR below x
                      if xR[i] > x+xR_max-domainLength_norm+xR[1] then
                         val_index = i-1                         
                         break
                      else 
                        -- If this is not true, then we have reached a point in x greater than the right profile
                        -- In this case, just set the value to the last value
                        val_index = numR
                      end
                   end
                   return n_e0R[val_index]*n_e 
                else 
                   return n_e
                end
              end,
              temperature = function (t, xn)
                 return T_e
              end,
              driftSpeed = function (t, xn)
                 return {vd_e}
              end
             },
      evolve = true, -- evolve species?
      bcx = { Plasma.AbsorbBC{},
               Plasma.AbsorbBC{} },
      diagnostics = { "M0", "M1i", "M2","M3i","intM0" },
      collE=Plasma.LBOCollisions {
         collideWith = {'elc','ion'},
         frequencies = {nu_eeProfile, nu_eiProfile},
      },
      src = Plasma.SteadySource {
         sourceSpecies = {"ion"},
         sourceLength = 1.0,
         profile = function (t, xn)
            local x, v = xn[1], xn[2]
               local m = maxwellian(1.0,0.0,vth_e,v)
               if x <= sheathLength then
                  return 0.0
               elseif x > sheathLength and x <= sheathLength + linearLength then
                  return fL(x)*m
               elseif x > sheathLength + linearLength and x <= sheathLength + linearLength + centerLength then
                  return m/(linearLength + centerLength)
               elseif x > sheathLength + linearLength + centerLength and x <= sheathLength + 2.0*linearLength + centerLength then
                  return fR(x)*m
               else
                  return 0.0
               end
            end,
      },
   },

   ion = Plasma.Species {
      charge = q_i, mass = m_i,
      -- velocity space grid
      lower = {-6.0*vth_i},
      upper = {6.0*vth_i},
      cells = {64},
      -- initial conditions
      init = Plasma.MaxwellianProjection {
              density = function (t, xn)
                -- get position normalized by Debye length 
                local x = xn[1]/lambda_D      
                
                -- Make an index for where we want to take the values from
                local val_index = 1
                
                -- Split this into three different regions, left, center, and right
                if x < xL_max then
                   -- We are in the left region                
                   -- Iterate through xL to determine where we should get the data from
                   for i = 1,numL do 
                      -- We want to use the value that corresponds to the highest value of xL below x
                      -- One error that comes up with this method is that some values of x initially might be negative
                      -- To fix this, we will set up an error catching early on to set and negative x value points to the first index
                      -- This is accomplished by skipping any x <= 0                                              
                      if xL[i] > x and x > 0 then
                         val_index = i-1                         
                         break
                      end
                   end
                   return n_i0L[val_index]*n_i
                elseif x > domainLength_norm-xR_max+xR[1] then
                  -- We are in the right region
                  -- Note that we had to add by xR[1] since the first value is not at 0 for the right profiles
                  for i = 1,numR do 
                      -- We want to use the value that corresponds to the highest value of xR below x
                      if xR[i] > x+xR_max-domainLength_norm+xR[1] then
                         val_index = i-1                         
                         break
                      else 
                        -- If this is not true, then we have reached a point in x greater than the right profile
                        -- In this case, just set the value to the last value
                        val_index = numR
                      end
                   end
                   return n_i0R[val_index]*n_i
                else 
                   return n_i
                end
              end,
              temperature = function (t, xn)
                 return T_i
              end,
              driftSpeed = function (t, xn)
                -- get position normalized by Debye length 
                local x = xn[1]/lambda_D      
                
                -- Make an index for where we want to take the values from
                local val_index = 1
                
                -- Split this into three different regions, left, center, and right
                if x < xL_max then
                   -- We are in the left region                
                   -- Iterate through xL to determine where we should get the data from
                   for i = 1,numL do 
                      -- We want to use the value that corresponds to the highest value of xL below x
                      -- One error that comes up with this method is that some values of x initially might be negative
                      -- To fix this, we will set up an error catching early on to set and negative x value points to the first index
                      -- This is accomplished by skipping any x <= 0                                              
                      if xL[i] > x and x > 0 then
                         val_index = i-1                         
                         break
                      end
                   end
                   return {uiL[val_index]*cs}
                elseif x > domainLength_norm-xR_max+xR[1] then
                  -- We are in the right region
                  -- Note that we had to add by xR[1] since the first value is not at 0 for the right profiles
                  for i = 1,numR do 
                      -- We want to use the value that corresponds to the highest value of xR below x
                      if xR[i] > x+xR_max-domainLength_norm+xR[1] then
                         val_index = i-1                         
                         break
                      else 
                        -- If this is not true, then we have reached a point in x greater than the right profile
                        -- In this case, just set the value to the last value
                        val_index = numR
                      end
                   end
                   return {uiR[val_index]*cs}
                else 
                   return {vd_i}
                end
              end
             },
      evolve = true, -- evolve species?
      bcx = { Plasma.AbsorbBC{},
               Plasma.AbsorbBC{} },
      diagnostics = { "M0", "M1i", "M2","M3i","intM0" },
      collI=Plasma.LBOCollisions {
         collideWith = {'elc','ion'},
         frequencies = {nu_ieProfile, nu_iiProfile},
      },
      src = Plasma.SteadySource {
         sourceSpecies = {"ion"},
         sourceLength = 1.0,
         profile = function (t, xn)
            local x, v = xn[1], xn[2]
               local m = maxwellian(1.0,0.0,vth_i,v)
               if x <= sheathLength then
                  return 0.0
               elseif x > sheathLength and x <= sheathLength + linearLength then
                  return fL(x)*m
               elseif x > sheathLength + linearLength and x <= sheathLength + linearLength + centerLength then
                  return m/(linearLength + centerLength)
               elseif x > sheathLength + linearLength + centerLength and x <= sheathLength + 2.0*linearLength + centerLength then
                  return fR(x)*m
               else
                  return 0.0
               end
            end,
      },
   },

   -- Field solver.
   field = Plasma.Field {
      epsilon0 = epsilon_0,
      evolve   = true, -- Evolve field?
      hasMagneticField = false,
      bcLowerPhi = {{ T ="D", V = V_left}},
      bcUpperPhi = {{ T ="D", V = 0.0}},
   },
}
-- run application
sim:run()


