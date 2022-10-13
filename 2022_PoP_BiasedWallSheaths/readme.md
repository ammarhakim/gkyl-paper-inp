2022 Continuum kinetic investigation of impact of bias potential on sheath formation in the FuZE parameter regime (PoP paper)

All results in this paper use revision [e53362c298c5](https://github.com/ammarhakim/gkyl/commit/e53362c298c5d57149c4d3c77e10b755fdcdb153).

1. All of the 1X-1V simulations are in the 1x1v_input_files directory. They span simulations with bias potentials from 0 to 10 kV. 
2. All of the 1X-2V simulations are in the 1x2v_input_files directory. The only input file is for the 5 kV case.
3. A python notebook, Robertson_solver.ipynb, is provided for setting initial conditions. The first cell is a function developed to numerically solve for the initial conditions based on [Robertson, S. (2013) PPCF](https://doi.org/10.1088/0741-3335/55/9/093001). The second cell saves the initial conditions for the left and right walls in text files. For this paper, the inputs for calculating the initial conditions are:
   phiL_SI: 0 to 10 kV
   Ti: 2000 eV
   Te: 2000 eV
   mi = 1.67262192369e-27 kg
   expected_niui = 0.55
   