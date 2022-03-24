-- Gkyl --------------`:`--------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Projection = require("App.Projection")

-- Constants
mass = 1.0
n0 = 0.5
grav = 1.0
alph = 25.0

Atwood = 1.0/3.0
Kn = 0.1

Lx = 0.75
Ly = 1.5

T0 = 2*grav*mass*(3*Ly*alph+Ly*math.log(math.cosh(alph)))/(3*alph)

-- Velocity space extents
T_center = 2*T0
vth_center = math.sqrt(T_center/mass)

epsilon0, mu0 = 1.0, 1.0

numCells = 64
kNumber = math.pi/(2*Lx)
tauRT = math.sqrt(1/(kNumber*Atwood*grav))

maxwellian = function(vx, vy, vz, n, ux, uy, uz, vth)
   return n/((math.sqrt(2*math.pi*vth^2))^3)*math.exp(-((vx-ux)^2 + (vy-uy)^2 + (vz-uz)^2)/(2*vth^2))
end

dy = Ly/numCells

diags = {"M0","M2ij","M3i","VtSq","Udrift"}

App = Vlasov.App{
    logToFile = false,

    tEnd = 3.0*tauRT,
    nFrame = 10,
    lower = {-Lx,-Ly},
    upper = {Lx, Ly},
    cells = {numCells,numCells},
    basis = "serendipity",
    polyOrder = 2,
    timeStepper = "rk3s4",
    decompCuts = {16,32},
    useShared = false,
    writeGhost = false,
    periodicDirs = {1},
    
    neut = Vlasov.Species {
        charge = 0.0, mass = 1.0,
        
        -- Velocity space grid
        lower = {-4.0*vth_center, -4.0*vth_center, -4.0*vth_center},
        upper = {4.0*vth_center, 4.0*vth_center, 4.0*vth_center},
        cells = {16,16,16},
     
		init = Vlasov.MaxwellianProjection {
			density = function (t, xn)
				local x = xn[1]
				local y = xn[2]
				local numDens = n0/2*math.tanh(alph*y/Ly) + 3/2*n0
				return numDens	
			end,
			driftSpeed = function (t, xn)
				local x, y = xn[1], xn[2]
				local k = kNumber
				local yr = Ly/10
				local perturb = -0.1*vth_center*math.cos(k*x)*math.exp(-(y*y)/(2*yr*yr))
				return {0.0, perturb, 0.0}
			end, 
			temperature = function (t, xn)
				local x =  xn[1]
				local y = xn[2]
				local numDens = n0/2*math.tanh(alph*y/Ly) + 3/2*n0
				local Press = -grav*mass*(n0/2*math.log(math.cosh(alph*y/Ly))*Ly/alph + 3/2*n0*y) + 3/2*n0*T0
				local temp = Press/(numDens)
				return temp
			end,
			isInit = true,
		},			
 
        evolve = true,
        diagnostics = diags,
		bcy = {
		function (t, z)
			local x = z[1]
			local y = z[2] - Ly - dy/2
			local numDens = n0/2*math.tanh(alph*y/Ly) + 3/2*n0
			local Press = -grav*mass*(n0/2*math.log(math.cosh(alph*y/Ly))*Ly/alph + 3/2*n0*y) + 3/2*n0*T0
			local temp = Press/(numDens)
			local vth = math.sqrt(temp/mass)
			return maxwellian(z[3], z[4], z[5], numDens, 0, 0, 0, vth)
		end,
		function (t, z)	
			--local x, y = z[1], z[2]
			local x = z[1]
			local y = z[2] + Ly + dy/2
			local numDens = n0/2*math.tanh(alph*y/Ly) + 3/2*n0
			local Press = -grav*mass*(n0/2*math.log(math.cosh(alph*y/Ly))*Ly/alph + 3/2*n0*y) + 3/2*n0*T0
			local temp = Press/(numDens)
			local vth = math.sqrt(temp/mass)
			return maxwellian(z[3], z[4], z[5], numDens, 0, 0, 0, vth)
		end,
		},
	evolveFnBC = false,
	vlasovExtForceFunc = function(t, xn)
		return 0.0, -grav, 0.0
	end,
        
	coll = Vlasov.BGKCollisions {
		collideWith = {"neut"},
		frequencies = {vth_center/(Lx*Kn)},
	},
    },
}

App:run()