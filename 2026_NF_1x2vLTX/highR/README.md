High recycling simulations
--------------------------

#[ Reference simulation
- gt0

#[ Resolution scan
- gt1: 32x16x12
- gt0: 64x16x12
- gt2: 128x16x12
- gt3: 64x8x12
- gt4: 64x32x12
- gt5: 64x16x6
- gt6: 64x16x24

#[ vmax scan (using gt0 as base)
- gt7  : 1.5 * vparmax, 64x24x12
- gt8  : 2.0 * vparmax, 64x32x12
- gt9  : (2/3) * vperpmax, 64x16x8
- gt10 : 1.5 * vperpmax, 64x16x18

#[ kperp scan (using gt0 as base)
- gt11: 0.025
- gt12: 0.05
- gt13: 0.1
- gt0 : 0.15
- gt14: 0.2
- gt15: 0.25
- gt16: 0.3
- gt17: 0.35

#[ Scan nuFrac (using gt0 as base)
- gt18: 0.0125
- gt19: 0.025
- gt20: 0.05
- gt21: 0.1
- gt22: 0.2
- gt23: 0.35
- gt24: 0.5
- gt0 : 1.0
- gt25: 2.0
- gt26: 2.5
- gt27: 3.0
- gt28: 6.0
- gt29: 12.0

#[ Boltzmann electrons
- gt30

#[ Constant B
- gt31
