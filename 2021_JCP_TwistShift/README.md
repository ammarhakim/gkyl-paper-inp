2021 DG interpolation and sheared (twist-shift) BCs JCP paper

Input files used to generate the results presented in section 4.

1. tests-2d:
  - twistShift-2x-yGaussian: shifts a y-Gaussian forward and backward.
  - twistShift-2x-yStep: shifts a rectangular function in y by half a cell
    length to illustrate diffusion.
  - twistShift-2x-ySineDiffusion: shift a y-sine forward and backward
    iteratively to measure diffusion.
2. tests-3d:
  - twistShift-3x-staticGaussian: static application of twist-shift BCs to a 3D
    field.
  - passiveAdvection-3x: passively advect a rectangular function in z with
    twist-shift BCs.
3. tests-5d:
  - twistShift-3x2v: static application of twist-shift BCs in z to a 5D field.
  - n10-es-adiabatic-global-deltaf-twist: linear ITG Cyclone electrostatic
    gyrokinetic simulation with the delta-f model in a global domain with twist-shift BCs.
  - n10-es-adiabatic-global-deltaf: the same as
    n10-es-adiabatic-global-deltaf-twist but with periodic BCs.
  - n10-es-adiabatic-global-deltaf-2pi: the same as
    n10-es-adiabatic-global-deltaf but with Lz=6*pi.
