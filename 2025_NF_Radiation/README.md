Input files for "A Kinetic Line-Driven Radiation Operator and Its Application to Gyrokinetics" by Jonathan Roeltgen et. al.

Resolution scan (section 4.1)
- resolution_scan/
  - Ar1_resscan_12x6.c: Uniform with 12 cells in v<sub>||</sub> and 6 cells in &mu;
  - Ar1_resscan_12x6.c: Uniform with 12 cells in v<sub>||</sub> and 6 cells in &mu;
  - Ar1_resscan_16x8.c: Uniform with 16 cells in v<sub>||</sub> and 8 cells in &mu;
  - Ar1_resscan_32x16.c: Uniform with 32 cells in v<sub>||</sub> and 16 cells in &mu;
  - Ar1_resscan_64x32.c: Uniform with 64 cells in v<sub>||</sub> and 32 cells in &mu;
  - Ar1_resscan_128x64.c: Uniform with 128 cells in v<sub>||</sub> and 64 cells in &mu;
  - Ar1_resscan_nonuniform_16x8.c: Non-uniform with 16 cells in v<sub>||</sub> and 8 cells in &mu;

Coronal Equilibrium (section 4.2)
- coronal/
  - C_2eV_coronal.c: All carbon charge states at 2eV
  - C_5eV_coronal.c: All carbon charge states at 5eV
  - C_7eV_coronal.c: All carbon charge states at 7eV
  - C_10eV_coronal.c: All carbon charge states at 10eV
  - C_20eV_coronal.c: All carbon charge states at 20eV
  - C_30eV_coronal.c: All carbon charge states at 30eV
  - C_50eV_coronal.c: All carbon charge states at 50eV
  - C_70eV_coronal.c: All carbon charge states at 70eV
  - C_100eV_coronal.c: All carbon charge states at 100eV
  - C_150eV_coronal.c: All carbon charge states at 150eV
  - C_200eV_coronal.c: All carbon charge states at 200eV
  - C_250eV_coronal.c: All carbon charge states at 250eV
  - C_300eV_coronal.c: All carbon charge states at 300eV

Non-Maxwellian comparison (section 5)
- SIKE_compare/
  - example_dist/ Contains the electron distributions, velocity grids, density and temperatures used. Originally calculated from SOL-KiT.
  - non_maxwellian_compare.c: reads in distribution functions from example_dist/ and uses line temperature from line specified with -w option to input file
  - run_all.sh: An example shell script to run non_maxwellian_compare.c for each available temperature
  - git_hash: The git hash it was run on.
