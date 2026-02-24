local Vlasov = G0.Vlasov
pi = math.pi  

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.   
T = 1.0 -- Temperature (units of mc^2)
injRate = 0.5 -- Injection rate.  

-- Simulation parameters.
Nx =  200 -- Cell count (configuration space: x-direction).
Nvx = 800 -- Cell count (velocity space: vx-direction).
Lx = 20.0 -- Domain size (configuration space: x-direction).
kx = 2.0*pi/Lx -- Smallest wavenumber (largest wavelength) in the domain.
vx_max = 10000.0  -- Domain boundary (electron velocity space: vx-direction).
alpha = 0.4 --1.0 -- 10.0 -- Applied perturbation amplitude.
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.2 -- CFL coefficient.

t_end = 1800.0 -- Final simulation time.
num_frames = 900 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (only x).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID  = G0.Model.SR,
    writeCellAvg = true,
    charge = charge_elc, mass = mass_elc,
    skipCellThresh = 1.0e-12,
    
    -- Velocity space grid.
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Nvx },

    mapc2pVel = {
      -- vx mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0

          local ncells = Nvx
          local v1, v2
          local pmin, pint, pmax, p0, vZero

          pmin = 0.14
          pmax = vx_max

          if (vc < 0.0) then
            v1 = pmin*vc*ncells-math.exp((-vc)*math.log(vx_max))+1
            vp = v1;
          else
            v1 = pmin*vc*ncells+math.exp(vc*math.log(vx_max))-1
            vp = v1;
          end
          return vp
        end
      },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
	 local x = xn[1]
	 math.randomseed(0)
	 local n = 1.0
	 for i = 1, 5 do
	   n = n + 0.5*alpha*math.random()*math.cos(kx*i*x  + 2.0 * pi * math.random())
	 end
	 if (n<0.0) then
	   n=0.0;
	 end	 
	 --if (n<0.0) then
	 --  n=0.0
	 --end
	 -- n = n + 0.1
	 return n  -- Electron total number density.
        end,
        temperatureInit = function (t, xn)
          return T -- Electron isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0 -- Electron drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.Proj,
      sourceLength = Lx,
      sourceSpecies = "elc",
 
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,
 
          densityInit = function (t, xn)
            local n = injRate
            return n
          end, 
 
 	  temperatureInit = function (t, xn)
            return T -- Electron source isotropic temperature.
          end,
 	  
          driftVelocityInit = function (t, xn)
            return 0.0 -- Electron source drift velocity.
          end
        }
      }
   },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_ion, mass = mass_ion,
    writeCellAvg = true,
    skipCellThresh = 1.0e-12,  
    
    -- Velocity space grid.
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Nvx },

    mapc2pVel = {
      -- vx mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0

          local ncells = Nvx
          local v1, v2
          local pmin, pint, pmax, p0, vZero

          pmin = 0.14
          pmax = vx_max

          if (vc < 0.0) then
            v1 = pmin*vc*ncells-math.exp((-vc)*math.log(vx_max))+1
            vp = v1;
          else
            v1 = pmin*vc*ncells+math.exp(vc*math.log(vx_max))-1
            vp = v1;
          end
          return vp
        end
      },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
	 local x = xn[1]
	 math.randomseed(0)
	 n = 1.0
	 for i = 1, 5 do
	   n = n - 0.5 * alpha*math.random()*math.cos(kx*i*x + 2.0 * pi * math.random())
	 end
	 if (n<0.0) then
	   n=0.0;
	 end
	 -- else
	 --  n = -n;
	 -- end
	 -- n = n+0.1
	 return n -- Ion total number density.
	 
        end,
        temperatureInit = function (t, xn)
          return T -- Ion isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0 -- Ion drift velocity.
        end
      }
    },

   source = {
      sourceID = G0.Source.Proj,
      sourceLength = Lx,
      sourceSpecies = "ion",
   
      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,
   
          densityInit = function (t, xn)
            return injRate
          end,
          temperatureInit = function (t, xn)
            return T -- Ion source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return 0.0 -- Ion source drift velocity.
          end
        }
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }

  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      Exnew = 400.0 -- Total electric field (x-direction).

      math.randomseed(0)  
      for i = 1, 5 do
         Exnew = Exnew - alpha*math.random()*math.sin(kx*i*x + 2.0 * pi * math.random())/(i*kx)
      end

      local Ex = Exnew
      
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }  
}

vlasovApp:run()