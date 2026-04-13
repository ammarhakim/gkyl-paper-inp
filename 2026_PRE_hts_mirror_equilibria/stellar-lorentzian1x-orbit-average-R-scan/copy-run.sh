#!/bin/bash

# Perlmutter-specific modules (commented out for stellar-amd)
# module load PrgEnv-gnu/8.5.0
# module load craype-accel-nvidia80
# module load cray-mpich/8.1.28
# module load cudatoolkit/12.4
# module load nccl/2.18.3-cu12

# Define arrays for magnetic field parameters
mcB_values=(2.130115 2.665626 3.691260 4.490901 5.416264)
gamma_values=(0.451454 0.331696 0.226381 0.182792 0.149893)
R_values=(3 5 10 15 22)
src_amp=(30657257.3675 30854916.3067 30975590.4361 30958371.5349 31019560.5328)
src_temp=(18907.9098163 18497.8733719 18282.7074073 18204.6444904 18149.4757841)

mkdir -p R-scan

rm -f core/sim

# Submit jobs for paired mcB and gamma scans
for i in "${!mcB_values[@]}"; do
  mcB="${mcB_values[$i]}"
  gamma="${gamma_values[$i]}"
  R="${R_values[$i]}"
  
  # Create the folder structure
  folder_name="R-scan/R-${R}"

  # mkdir -p "$folder_name"

  # Copy core files into the folder
  cp core/* "$folder_name/"

  # Change into the folder
  cd "$folder_name" || exit

	if [ "$R" -eq 3 ]; then
		sed -i "152s/.*/  int Nz = 192;/" sim.c
	elif [ "$R" -eq 5 ]; then
		sed -i "152s/.*/  int Nz = 256;/" sim.c
	fi

  sed -i "160s/.*/  double ion_source_amplitude = ${src_amp[$i]};/" sim.c
  sed -i "161s/.*/  double ion_source_temp = ${src_temp[$i]} * eV;/" sim.c
  sed -i "167s/.*/  double mcB = $mcB;/" sim.c
  sed -i "168s/.*/  double gamma = $gamma;/" sim.c
  sed -i "478s|.*|    .filename_psi = \"/home/mr1884/scratch/gkylmax/generate_efit/lorentzian_R${R}.geqdsk_psi.gkyl\",|" sim.c
	
  sed -i "5s/.*/#SBATCH -J poa-R-${R}/" jobscript-gkyl-stellar-amd
  
  # Build the simulation
	make clean
  make

  # ./sim -s1

  # Submit the job
  # bash optimize_source_params.sh
  sbatch jobscript-gkyl-stellar-amd
	# bash submit-restarts.sh
	wait

  # Print confirmation
  echo "submitted job for R = $R (mcB = $mcB, gamma = $gamma)"

  # Change back to the root directory
  cd - || exit
done