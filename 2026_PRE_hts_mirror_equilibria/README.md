# 2026_XX_high_field_mirror_equilibria

These files are for the publication:
M. Rosen, M. Francisquez, A. Hakim, G. W. Hammett. "Gyrokinetic equilibria of high temperature superconducting magnetic mirrors" (2026).

Simulations are run on commit https://github.com/ammarhakim/gkeyll/commit/eacf83eeea0c7eef2e8d2e68c9c548665c7c2262

Note, you need to modify the shared makefile to include the common header folder.

## What is in this directory?

This folder is organized around four tasks:

1. Run 1x-2v gyrokinetic simulations (single run and parameter scans).
2. Helper folders (common-header, generate_efit)
3. Post-process outputs and compare with analytical confinement models.

## Directory map

- `analysis/`
	- Post-processing and validation scripts (e.g., confinement scaling, distribution-function diagnostics, moment comparisons).
	- Representative scripts:
		- `confinement-time.py`: mirror-ratio scaling fits for Maxwellian/beam source cases.
		- `calc_pastukhov_potential.py`: infer $e\phi/T_e$ from simulation-derived confinement times.
		- `validation.py`: POA vs FDP profile comparison plots.
		- `distribution_functions.py`: velocity-space diagnostics at selected axial locations.

- `common-header/`
	- Shared C header + build/job helpers used across simulation directories.
	- `sim.h` defines the common POA phase state/parameter structs and utility routines.

- `generate_efit/`
	- Geometry generation for double-Lorentzian mirrors.
	- `write_efit_double_lorentzian.py` writes `lorentzian_R*.geqdsk` files for selected mirror ratios.
	- `*.geqdsk_psi.gkyl` files are pre-generated geometry inputs consumed by simulations.

- `stellar-lorentzian1x-kinetic-exploration/`
	- Produces the kinetic electron simulation. Takes input as the output of `stellar-lorentzian1x-orbit-average-beams` (`sim.c`, `Makefile`, job script).

- `stellar-lorentzian1x-orbit-average/`
	- Main  Maxwellian-source run setup.
	- Contains `sim.c`, build scripts, plotting helper, and restart submission scripts.

- `stellar-lorentzian1x-orbit-average-validate/`
	- Validation run configuration used for POA/FDP cross-checks.

- `stellar-lorentzian1x-orbit-average-beams/`
	- Beam-source variant of the orbit-average setup.
	- Includes `plot.sh` for quick `pgkyl` diagnostics and figure generation.

- `stellar-lorentzian1x-orbit-average-beams-validate/`
	- Beam-source validation cases and comparison plotting scripts.

- `stellar-lorentzian1x-orbit-average-R-scan/`
	- Maxwellian-source mirror-ratio scan infrastructure.
	- `core/` stores base templates.
	- `copy-run.sh` stamps out `R-scan/R-*` cases and updates per-case `sim.c` parameters.
	- `find-geo-coeffs.py` solves for Lorentzian coefficients at target mirror ratios.

- `stellar-lorentzian1x-orbit-average-beams-R-scan/`
	- Beam-source mirror-ratio scan infrastructure.
	- Includes tools for source tuning (`optimize_source_params.py`) and confinement summaries (`info.sh`, `pastukhov_calc.py`).


## Typical workflow

1. Build Gkeyll, checkout this specific commit, modify shared makefile to include common header
2. Build simulation executable in a run directory.
	 - Example: `cd stellar-lorentzian1x-orbit-average && make`
3. Launch the run (directly or with scheduler scripts).
	 - Example cluster flow: `sbatch jobscript_gkyl-stellar-amd`
4. For R scans, create/submit case directories from `core/` templates.
	 - Example: run `copy-run.sh` inside an `*-R-scan/` directory.
5. Post-process with `analysis/*.py` or local `plot.sh` helpers.