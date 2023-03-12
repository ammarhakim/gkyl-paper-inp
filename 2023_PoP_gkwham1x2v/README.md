2023 PoP paper on 1x2v gyrokinetic simulations of WHAM

Input files used to generate results in the paper about 1x2v gyrokinetic simulations of the Wisconsin HTS Axisymmetric Mirror (WHAM). 

1. Boltzmann electron simulations.
  - gk57-wham1x2v: base simulation with 244x64x192 resolution and no force-softening.
  - gk71-wham1x2v: simulation with 244x64x192 resolution and force-softening.
  - gk77-wham1x2v: simulation with 244x96x192 resolution and force-softening.
  - resolution scan:
    * gk57-wham1x2v: 244x64x192
    * gk60-wham1x2v: 244x128x192
    * gk63-wham1x2v: 244x256x192
    * gk61-wham1x2v: 244x64x320
2. Kinetic electron simulations with linear polarization.
  - kperp * rhos scan:
    * gk52-wham1x2v: 0.005
    * gk48-wham1x2v: 0.01
    * gk46-wham1x2v: 0.05
    * gk41-wham1x2v: 0.1
    * gk45-wham1x2v: 0.15
    * gk44-wham1x2v: 0.2
    * gk43-wham1x2v: 0.25
    * gk32-wham1x2v: 0.3
  - resolution scan:
    * gk54-wham1x2v: 288x64x32
    * gk49-wham1x2v: 288x96x32
    * gk48-wham1x2v: 288x192x32
    * gk50-wham1x2v: 288x384x32
    * gk51-wham1x2v: 288x192x64
    * gk53-wham1x2v: 288x192x128
    * gk55-wham1x2v: 288x64x192
  - gk55-wham1x2v: base simulation
3. Kinetic electron simulations with nonlinear polarization and force-softening.
  - gk75-wham1x2v 
