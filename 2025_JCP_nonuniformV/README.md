2025 Nonuniform velocity discretization paper

Input files used to generate the results presented in section 4 of the paper titled "Conservative velocity mappings for discontinuous Galerkin kinetics".

1. Collisionless damping of an ion acoustic wave:
  - gk_ion_sound_adiabatic_elc_uni_1x2v_p1.c: k rho=0.125 test with a uniform velocity grid.
  - gk_ion_sound_adiabatic_elc_nonuni_1x2v_p1.c: k rho=0.125 test with a nonuniform velocity grid.
2. Collisional isotropization:
  - gk_lbo_relax_uni_1x2v_p1.c: test with a uniform velocity grid.
  - gk_lbo_relax_nonuni_1x2v_p1.c: test with a nonuniform velocity grid.
3. 1D HTS mirror:
  - gz57_1x2v_p1.c: 1D simulation of a high-temperature superconducting mirror
    with nonuniform velocity grid.
4. 2D ASDEX Upgrade:
  - gk_aug36190_sol_2x2v_p1.c: 2D simulation of the ASDEX Upgrade SOL with
    nonuniform velocity grid.
3. 3D LAPD:
  - gk_lapd_nonuniv_3x2v_p1.c: 3D simulation of LAPD turbulence with nonuniform
    velocity grid.
