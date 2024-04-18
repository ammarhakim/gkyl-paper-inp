-- Gkyl --------------`:`--------------------------------------------------------
local Vlasov = require("App.PlasmaOnCartGrid").VlasovMaxwell()
local Projection = require("App.Projection")

-- Constants
local mass = 1.0
local n0 = 0.5
local grav = 1.0
local alph = 25.0

local Atwood = 1.0 / 3.0

-- Kn defined here to set normalized coll freq
local Kn = 0.01

local fac = 1.0
local Lx = 0.5
local Ly = 1.0

local T0 = fac * 2 * grav * mass * (3 * Ly * alph + Ly * math.log(math.cosh(alph))) / (3 * alph)

-- Velocity space extents
local T_center = 2 * T0
local vth_center = math.sqrt(T_center / mass)
-- local normNu = 500.0

local epsilon0, mu0 = 1.0, 1.0

local numCells = 32
local kNumber = math.pi / (2 * Lx)
local tauRT = math.sqrt(1 / (kNumber * Atwood * grav))

maxwellian = function(vx, vy, vz, n, ux, uy, uz, vth)
    return n / ((math.sqrt(2 * math.pi * vth ^ 2)) ^ 3) *
               math.exp(-((vx - ux) ^ 2 + (vy - uy) ^ 2 + (vz - uz) ^ 2) / (2 * vth ^ 2))
end

dy = (2 * Ly) / numCells
dx = (2 * Lx) / numCells

-- diags = {"M0", "M2ij"}
-- -- print(tauRT)
-- App = Vlasov.App {
--     logToFile = true,

--     -- tEnd = 0.002413,
--     tEnd = 3.0 * tauRT,
--     nFrame = 15,
--     lower = {-Lx, -Ly},
--     upper = {Lx, Ly},
--     cells = {numCells, 2 * numCells},
--     basis = "serendipity",
--     polyOrder = 2,
--     timeStepper = "rk3s4",
--     decompCuts = {32, 64},
--     useShared = false,
--     writeGhost = false,
--     periodicDirs = {1},

--     neut = Vlasov.Species {
--         charge = 0.0,
--         mass = 1.0,

--         -- Velocity space grid
--         lower = {-5.0 * vth_center, -5.0 * vth_center, -5.0 * vth_center},
--         upper = {5.0 * vth_center, 5.0 * vth_center, 5.0 * vth_center},
--         cells = {16, 16, 16},

--         -- init = Projection.KineticProjection.FunctionProjection {
--         --     	func = function (t, xn) return end,
--         --     	fromFile = "neut_10.bp",
--         -- },
--         init = Vlasov.MaxwellianProjection {
--             density = function(t, xn)
--                 local x = xn[1]
--                 local y = xn[2]
--                 local numDens = n0 / 2 * math.tanh(alph * y / Ly) + 3 / 2 * n0
--                 return numDens
--             end,
--             driftSpeed = function(t, xn)
--                 local x, y = xn[1], xn[2]
--                 local k = kNumber
--                 local yr = Ly / 10
--                 local perturb = -0.1 * vth_center * math.cos(k * x) * math.exp(-(y * y) / (2 * yr * yr))
--                 return {0.0, perturb, 0.0}
--             end,
--             temperature = function(t, xn)
--                 local x = xn[1]
--                 local y = xn[2]
--                 local numDens = n0 / 2 * math.tanh(alph * y / Ly) + 3 / 2 * n0
--                 local Press =
--                     -grav * mass * (n0 / 2 * math.log(math.cosh(alph * y / Ly)) * Ly / alph + 3 / 2 * n0 * y) + 3 / 2 *
--                         n0 * T0
--                 local temp = Press / (numDens)
--                 return temp
--             end,
--             isInit = true
--         },

--         vlasovExtForceFunc = function(t, xn)
--             return 0.0, -grav, 0.0
--         end,

--         evolve = true,
--         diagnostics = diags,
--         -- bcy = {Vlasov.ReflectBC{}, Vlasov.ReflectBC{}},
--         bcy = {function(t, z)
--             -- local x, y = z[1], z[2]
--             local x = z[1]
--             local y = z[2] - Ly - dy / 2
--             local numDens = n0 / 2 * math.tanh(alph * y / Ly) + 3 / 2 * n0
--             local Press =
--                 -grav * mass * (n0 / 2 * math.log(math.cosh(alph * y / Ly)) * Ly / alph + 3 / 2 * n0 * y) + 3 / 2 * n0 *
--                     T0
--             local temp = Press / (numDens)
--             local vth = math.sqrt(temp / mass)
--             return maxwellian(z[3], z[4], z[5], numDens, 0, 0, 0, vth)
--         end, function(t, z)
--             -- local x, y = z[1], z[2]
--             local x = z[1]
--             local y = z[2] + Ly + dy / 2
--             local numDens = n0 / 2 * math.tanh(alph * y / Ly) + 3 / 2 * n0
--             local Press =
--                 -grav * mass * (n0 / 2 * math.log(math.cosh(alph * y / Ly)) * Ly / alph + 3 / 2 * n0 * y) + 3 / 2 * n0 *
--                     T0
--             local temp = Press / (numDens)
--             local vth = math.sqrt(temp / mass)
--             return maxwellian(z[3], z[4], z[5], numDens, 0, 0, 0, vth)
--         end},
--         evolveFnBC = false,

--         coll = Vlasov.BGKCollisions {
--             collideWith = {"neut"},
--             -- frequencies = {vth_center/(Lx*Kn)},
--             nuFrac = 2500.0,
--             normNu = {1.0}
--         }
--     }

-- }

-- App:run()

require "Diags"

for frame = 0, 15 do
    local tbl = {
        lower = {-Lx, -Ly},
        upper = {Lx, Ly},
        vExt = 5.0 * vth_center,
        numCells = {32, 64},
        decompCuts = {32, 32},
        frame = frame,
        -- normalize = true,
        -- factor = n0 * vth_center ^ 3,
        writeVtSq = true
    }

    local tbl = InitTbl(tbl)
    calcEnergyFluxTerms(tbl)
    calcNonMaxDensity(tbl)

end
