-- Gkyl ------------------------------------------------------------------------
--
-- Shift a y-Gaussian profile, then shift it back. Do this a number of times
-- to evaluate how diffusive the algorithm is.
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
local numCells        = {1, 160}
local periodicDirs    = {2}
local yShiftPolyOrder = polyOrder

local Ly  = upper[2]-lower[2]
local dy  = Ly/numCells[2]
local ky  = (2*math.pi)/Ly

-- Number of iterations. Each iteration consists of
-- one forward and one backward shift.
local numIter = 1000

local nFrames = 20

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

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
                      return 0.9+2.*dy/3.
                   end

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function(t, xn) return 1. end,  -- Set later.
}
-- Donor field function.
local fldDoFunc = function(t, xn)
   local x, y = xn[1], xn[2]
   return 2.+math.cos(ky*y)
end

function fileName(varNm, fr)
  return string.format(varNm .. "_Nx%dNy%dP%d_yShP%deq1p1_%d.bp",grid:numCells(1),grid:numCells(2),polyOrder,yShiftPolyOrder,fr)
end

-- Project donor field function onto basis.
project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
project:advance(0., {}, {fldDo})

local cFr = 0
fldDo:write(fileName("fldDo", cFr), 0., cFr)
fldTar:write(fileName("fldTar", cFr), 0., cFr)

local intQuant = Updater.CartFieldIntegratedQuantCalc {
   onGrid = grid,   numComponents = 1,
   basis  = basis,  quantity      = "V",
}
local intFldDo  = DataStruct.DynVector { numComponents = 1, }
local intFldTar = DataStruct.DynVector { numComponents = 1, }

-- Updater for forward shift.
local twistShiftUpd = Updater.TwistShiftBC {
   onGrid = grid,   yShiftFunc      = yShiftFunc, 
   basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
}

-- Updater for backward shift.
local yShiftBackFunc = function(t, xn) return -yShiftFunc(t, xn) end
local twistShiftBackUpd = Updater.TwistShiftBC {
   onGrid = grid,   yShiftFunc      = yShiftBackFunc, 
   basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
}

local t1 = os.clock()
local shiftT2 = 0.
for it = 1, numIter do
   fldTar:clear(0.)
   local shiftT1 = os.clock()
   twistShiftUpd:advance(0., {fldDo}, {fldTar})   -- Forward shift.
   shiftT2 = shiftT2 + os.clock() - shiftT1

--   if (it % (numIter/nFrames) == 0) then fldTar:write(fileName("fldTar", it), it, it) end

   fldDo:clear(0.)
   local shiftT1 = os.clock()
   twistShiftBackUpd:advance(0., {fldTar}, {fldDo})
   shiftT2 = shiftT2 + os.clock() - shiftT1

   if (it % (numIter/nFrames) == 0) then
      cFr = cFr+1
      fldDo:write(fileName("fldDo", cFr), it, cFr)
   end

   intQuant:advance(it, {fldDo}, {intFldDo})
--   intQuant:advance(it, {fldTar}, {intFldTar})
end
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")
io.write("Time per shift: ", shiftT2/(2*numIter), " s\n")

intFldDo:write(fileName("intFldDo",numIter), numIter, numIter)
intFldTar:write(fileName("intFldTar",numIter), numIter, numIter)

