Input files used to simulate AUG discharge \#36190 (adapted from https://arxiv.org/pdf/2505.10754v3), and ITER SOLPS case \#2292, with a static neutral argon background added to study impurity transport.

The simulations were performed in branch `diffusive_bcs_mod` in commit `c8ded37fc`
# Base simulations

### these simulations are run from the output of the previous

ITER:

- ITER_no_impurities.c: The main-plasma simulation, without impurity argon added.
- ITER_impurities_static_ion.c: The simulation with impurity argon added, where the electrons and main ions are kept static to speed up convergence. NOTE: this script was lost, and attempted to be recreated post facto. Correctness of this script cannot be guaranteed.
- ITER_impurities_static_elc.c: Same as above, except now only the electrons are kept static, and the ions are evolved.
- ITER_impurities_dynamic_elc.c: The full simulation, with the electrons evolved as well (all physics enabled).

ASDEX:

- ASDEX_no_impurities.c: The main-plasma simulation, without impurity argon added.
- ASDEX_impurities_static_elc.c: The simulation with impurity argon added, where the electrons are kept static to speed up convergence
- ASDEX_impurities_dynamic_elc.c: The full simulation, with the electrons evolved as well (all physics enabled).

# Diffusivity scan

### A scan of the diffusivity of the impurity species. The diffusivity of the main plasma species are kept constant

- The number at the end indicates the diffusivity (i.e. 005 => D=0.05, 16 => D=1.6, etc.)

# Impurity species scan

### A variation of the impurity species from Ar to N and O respectively.

# Final remarks

- To my knowledge, the baseline EQDSK for ITER is not publically available, and therefore not included here, but can be obtained from the ITER baseline sharepoint (SOLPS_EQDSK.txt).

- All Gkeyll output files are available on request at the Dutch Institute For Fundamental Energy Research (DIFFER)
