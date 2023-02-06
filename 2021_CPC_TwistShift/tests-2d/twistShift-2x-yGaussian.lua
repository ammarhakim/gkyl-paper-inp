-- Gkyl ------------------------------------------------------------------------
--
-- Test the interpolation needed for twist-shift BCs in gyrokinetics in 2D.
--
-- Create two fields on a 2D grid, the donor field and the target field. Then
-- shift the donor field with a shift in y (that may be a function of x)
-- assuming periodicity in y to obtain the target field.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local polyOrder       = 1 
local lower           = {-2.0, -1.50}
local upper           = { 2.0,  1.50}
local numCells        = {1, 10}
local periodicDirs    = {2}
local yShiftPolyOrder = polyOrder

-- Parameters of the y-Gaussian.
local gauss_mu_y  = 0.
local gauss_sig_y = 0.3

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
                      return 1.1
                   end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
   }
   return fld
end

local wrapNum = function (val, lims, pickUpper)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   local lower, upper = lims.lo, lims.up
   local L        = upper - lower
   local disp     = (val - lower) % L
   local newCoord = lower + (L + disp) % L
   local eps      = 1.e-12
   if ( (lower-eps < newCoord and newCoord < lower + eps) or
        (upper-eps < newCoord and newCoord < upper + eps) ) then
      if pickUpper then 
         return upper
      else
         return lower
      end
   else
      return newCoord
   end
end

local grid = Grid.RectCart {
   lower = lower,  cells        = numCells,
   upper = upper,  periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }

local fldDo        = createField(grid, basis)
local fldDoShifted = createField(grid, basis)
local fldTar       = createField(grid, basis)

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn) return 1. end,  -- Set later.
}
-- Donor field function.
local fldDoFunc = function(t, xn)
   local x, y = xn[1], xn[2]
   local muY  = gauss_mu_y
   local sigY = gauss_sig_y
   return (1./math.sqrt(2.*math.pi*(sigY^2)))*math.exp(-((y-muY)^2)/(2.*(sigY^2)))
end
-- Shifted donor field function.
local fldDoShiftedFunc = function(t, xn)
   local x, y = xn[1], xn[2]
   local muY  = gauss_mu_y
   local sigY = gauss_sig_y
   return (1./math.sqrt(2.*math.pi*(sigY^2)))
      *math.exp(-((wrapNum(y-yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true)-muY)^2)/(2.*(sigY^2)))
end

function fileName(varNm)
  return string.format(varNm .. "_Nx%dNy%dP%d_yShP%deq1p1.bp",grid:numCells(1),grid:numCells(2),polyOrder,yShiftPolyOrder)
end

-- Project donor field function onto basis.
project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
project:advance(0., {}, {fldDo})
fldDo:write(fileName("fldDo"))
-- Project shifted donor field function onto basis.
project:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
project:advance(0., {}, {fldDoShifted})
fldDoShifted:write(fileName("fldDoShifted"))

local intQuant = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,   numComponents = 1,
   basis  = basis,  quantity      = "V",
}
local intFldDo = DataStruct.DynVector { numComponents = 1, }
intQuant:advance(0., {fldDo}, {intFldDo})
intFldDo:write(fileName("intFldDo"), 0., 0)

local twistShiftUpd = Updater.TwistShiftBC {
   onGrid = grid,   yShiftFunc      = yShiftFunc, 
   basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
}

local t1 = os.clock()
twistShiftUpd:_advance(0., {fldDo}, {fldTar})
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")

fldTar:write(fileName("fldTar"))

local intFldTar = DataStruct.DynVector { numComponents = 1, }
intQuant:advance(0., {fldTar}, {intFldTar})
intFldTar:write(fileName("intFldTar"), 0., 0)


-- ............... SHIFT BACK .................. --
local yShiftBackFunc = function(t, xn) return -yShiftFunc(t, xn) end

local twistShiftBackUpd = Updater.TwistShiftBC {
   onGrid = grid,   yShiftFunc      = yShiftBackFunc, 
   basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
}

fldDo:clear(0.)
local t1 = os.clock()
twistShiftBackUpd:_advance(0., {fldTar}, {fldDo})
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")

fldDo:write(fileName("fldDoBack"))

local intFldDoBack = DataStruct.DynVector { numComponents = 1, }
intQuant:advance(0., {fldDo}, {intFldDoBack})
intFldDoBack:write(fileName("intFldDoBack"), 0., 0)


