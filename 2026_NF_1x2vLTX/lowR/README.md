Low recycling simulations
--------------------------

#[ Reference simulation
- gr0

#[ Resolution scan
- gr1: 32x16x12
- gr0: 64x16x12
- gr2: 128x16x12
- gr3: 64x8x12
- gr4: 64x32x12
- gr5: 64x16x6
- gr6: 64x16x24

#[ vmax scan (using gr0 as base)
- gr7  : 1.5 * vparmax, 64x24x12
- gr8  : 2.0 * vparmax, 64x32x12
- gr9  : (2/3) * vperpmax, 64x16x8
- gr10 : 1.5 * vperpmax, 64x16x18

#[ kperp scan (using gr0 as base)
- gr11: 0.025
- gr12: 0.05
- gr13: 0.1
- gr0 : 0.15
- gr14: 0.2
- gr15: 0.25
- gr16: 0.3
- gr17: 0.35

#[ Scan nuFrac (using gr0 as base)
- gr18: 0.0125
- gr19: 0.025
- gr20: 0.05
- gr21: 0.1
- gr22: 0.2
- gr23: 0.35
- gr24: 0.5
- gr0 : 1.0
- gr25: 2.0
- gr26: 2.5
- gr27: 3.0
- gr28: 6.0
- gr29: 12.0

#[ Boltzmann electrons
- gr30

#[ Constant B
- gr31

#[ Constant B with Boltzmann electrons
- gr31_boltz_elc

#[ Like gr31 but higher resolution and larger vpar_max
- gr31_wo_transient_boltz_elc_vparMax2x_Nz2x_Nvpar8x

#[ Constant B, Boltzmann electrons, stronger sources to match gr0 steady state
- gr31_boltz_elc_match_gr0_Tsrc2x

#[ Scan collisionality in gr31_boltz_elc_match_gr0
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac0p0625: nu_frac=0.0625
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac0p125 : nu_frac=0.125
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac0p25  : nu_frac=0.25
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac0p5   : nu_frac=0.5
- gr31_boltz_elc_match_gr0_Tsrc2x             : nu_frac=1
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac2p0   : nu_frac=2
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac4p0   : nu_frac=4
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac8p0   : nu_frac=8
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac16p0  : nu_frac=16
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac32p0  : nu_frac=32
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac64p0  : nu_frac=64
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac128p0 : nu_frac=128
- gr31_boltz_elc_match_gr0_Tsrc2x_nuFrac256p0 : nu_frac=256

#[ Scan Te0 in gr31_boltz_elc_match_gr0 (Te0=178 eV)
- gr31_boltz_elc_match_gr0_Tsrc2x_Te0p25x: Te0 * 0.25
- gr31_boltz_elc_match_gr0_Tsrc2x_Te0p3x : Te0 * 0.3
- gr31_boltz_elc_match_gr0_Tsrc2x_Te0p4x : Te0 * 0.4
- gr31_boltz_elc_match_gr0_Tsrc2x_Te0p5x : Te0 * 0.5
- gr31_boltz_elc_match_gr0_Tsrc2x        : Te0 * 1
- gr31_boltz_elc_match_gr0_Tsrc2x_Te1p25x: Te0 * 1.25
- gr31_boltz_elc_match_gr0_Tsrc2x_Te1p5x : Te0 * 1.5
- gr31_boltz_elc_match_gr0_Tsrc2x_Te1p75x: Te0 * 1.75
- gr31_boltz_elc_match_gr0_Tsrc2x_Te2p0x : Te0 * 2
