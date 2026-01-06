2025 Nonuniform velocity discretization paper

Input files used to generate the results presented in section 4 of the paper titled "Conservative velocity mappings for discontinuous Galerkin kinetics". These simulations used [version eb83d1f of Gkeyll](https://github.com/ammarhakim/gkeyll/tree/eb83d1f6c0cf6e9c160ee637caaee8a2e30a9f64).

1. ion_acoustic_damping_linear: Linear collisionless damping of an ion acoustic wave.
  - uni: unfiorm velocity.
    * wavek0p125: k rho=0.125
    * wavek0p25 : k rho=0.25
    * wavek0p375: k rho=0.375
    * wavek0p5  : k rho=0.5
    * wavek0p75 : k rho=0.75
  - nonuni: nonuniform velocity.
    * wavek0p125: k rho=0.125
    * wavek0p25 : k rho=0.25
    * wavek0p375: k rho=0.375
    * wavek0p5  : k rho=0.5
    * wavek0p75 : k rho=0.75
  - nonuni_reduced_vparmax: nonuniform velocity, with reduced vparmax.
    * wavek0p125: k rho=0.125
    * wavek0p25 : k rho=0.25
    * wavek0p375: k rho=0.375
    * wavek0p5  : k rho=0.5
    * wavek0p75 : k rho=0.75
2. ion_acoustic_damping_nonlinear: Nonlinear collisionless damping of an ion acoustic wave.
  - isl0 : uniform vmax=6, 256x256x16
  - isl7 : 32x32x4
  - isl10: 32x128x4
  - isl6 : 64x64x4
  - isl9 : 64x128x4
  - isl15: 64x256x4
  - isl8 : 128x64x4
  - isl12: 128x128x4
  - isl14: 128x128x8
  - isl13: 128x256x4
  - isl3 : 256x256x32
  - isl2 : 256x512x16
  - isl1 : 512x256x16
  - isl4 : vpar_max=9, 384x256x16
  - isl5 : mu_max\*1.5, 256x384x16
  - isl11: k=0.5
  - isl16: isl13 vparmax=4, 128x172x4
  - isl17: vparmax=4.5, 128x192x4,
  - isl20: vperpmax=2.5, 128x256x4,
  - isl18: vperpmax=3, 128x256x4,
  - isl19: vperpmax=4.5, 128x256x4,
  - isl21: vparmax=4.5, vperpmax=2.5, 128x192x4
  - isl22: vparmax=5, vperpmax=5, 128x214x4
  - isl23: isl22, nonuniform vpar, 128x152x4
  - isl24: 128x64x4
  - isl25: 128x96x4
  - isl27: 128x108x4
  - isl26: 128x128x4
3. Collisional isotropization:
  - glb0 : uniform 1x12x16
  - glb1 : nonuni 1x4x4
  - glb2 : nonuni 1x6x4
  - glb3 : nonuni 1x12x4
  - glb4 : nonuni 1x18x4
  - glb5 : nonuni 1x24x4
  - glb6 : nonuni 1x4x6
  - glb7 : nonuni 1x6x6
  - glb8 : nonuni 1x12x6
  - glb9 : nonuni 1x18x6
  - glb10: nonuni 1x24x6
  - glb11: nonuni 1x4x12
  - glb12: nonuni 1x6x12
  - glb13: nonuni 1x12x12
  - glb14: nonuni 1x18x12
  - glb15: nonuni 1x24x12
  - glb16: nonuni 1x4x18
  - glb17: nonuni 1x6x18
  - glb18: nonuni 1x12x18
  - glb19: nonuni 1x18x18
  - glb20: nonuni 1x24x18
  - glb21: nonuni 1x4x24
  - glb22: nonuni 1x6x24
  - glb23: nonuni 1x12x24
  - glb24: nonuni 1x18x24
  - glb25: nonuni 1x24x24
  - glb26: nonuni 1x12x16
4. 1D HTS mirror:
  - gm0: nonuniform v
  - gm1: uniform v, same vmax as gm0
  - gm2: uniform v, vmax of linear part of gm0
5. 2D ASDEX Upgrade:
  - gaug11: D=0.2, 16x16x16x8
  - gaug12: D=0.3
  - gaug13: D=0.4
  - gaug14: gaug12, 32x16x16x8 = 56
  - gaug15: 64x16x16x8 = 112
  - gaug16: 16x32x16x8 = 56
  - gaug17: 32x32x16x8 = 112
  - gaug18: 64x32x16x8 = 224
  - gaug19: 16x64x16x8 = 448
  - gaug20: 32x64x16x8 = 896
  - gaug21: 64x64x16x8 = 1792
  - gaug22: 16x16x8x4  = 7
  - gaug23: 16x16x16x4 = 14
  - gaug24: 16x16x32x4 = 28
  - gaug25: 16x16x8x8  = 14
  - gaug26: 16x16x32x8 = 56
  - gaug27: 16x16x8x16 = 28
  - gaug28: 16x16x16x16= 56
  - gaug29: 16x16x32x16= 112
  - gaug30: uniform v
  - gaug31: uniform v, 16x16x66x64
  - gaug32: uniform v, vparmax of linear part of gaug12, 16x16x32x64
6. 3D LAPD:
  - gl0 : no bias, 36x36x8x8x4
  - gl1 : 24x24x8x8x4
  - gl2 : 72x72x8x8x4
  - gl3 : 128x128x8x8x4
  - gl4 : 36x36x16x8x4
  - gl5 : 36x36x24x8x4
  - gl6 : uniform v,
  - gl7 : uniform v, 36x36x8x12x16
  - gl8 : uniform v, 36x36x8x6x16, vmax of linear part of gl0
  - gl9 : bias=-7.275
  - gl10: bias=-3.075
  - gl11: bias=2.025
  - gl12: bias=9.825
