#!/bin/bash

# Define arrays for magnetic field parameters

mcB_values=(2.130115 2.665626 3.691260 4.490901 5.416264 6.51292)
gamma_values=(0.451454 0.331696 0.226381 0.182792 0.149893 0.124904)
R_values=(3 5 10 15 22 32)

# Initialize arrays to store optimized parameters
E_beam_midmap_values=()
T_beam_midmap_values=()
gamma_midmap_values=()

cp plot-comparison-bounce-average-vs-contour-integral.py temp.py

# Process each R value
for i in "${!R_values[@]}"; do
  R="${R_values[$i]}"
  mcB="${mcB_values[$i]}"
  gamma="${gamma_values[$i]}"
  
	sed -i "s/     mcB=.*/     mcB=$mcB,/" temp.py
	sed -i "s/   gamma_geo=.*/   gamma_geo=$gamma,/" temp.py

	output=$(python temp.py 2>&1)
	
	# Extract optimized parameters
	echo ""
	echo "=========================================="
	echo "EXTRACTED PARAMETERS for R = $R"
	echo "=========================================="
	
	# Extract E_beam_midmap (optimized) - field 5 is the optimized value after ->
	E_beam_opt=$(echo "$output" | grep "E_beam_midmap:.*->" | tail -1 | awk '{print $5}')
	
	# Extract T_beam_midmap (optimized) - field 5 is the optimized value after ->
	T_beam_opt=$(echo "$output" | grep "T_beam_midmap:.*->" | tail -1 | awk '{print $5}')
	
	# Extract gamma_midmap (optimized) - field 4 is the optimized value after ->
	gamma_midmap_opt=$(echo "$output" | grep "gamma_midmap:.*->" | tail -1 | awk '{print $4}')
	
	# Extract gamma_spatial (unchanged)
	gamma_spatial=$(echo "$output" | grep "gamma_spatial:" | tail -1 | awk '{print $2}')
	
	# Extract power comparisons
	power_spatial=$(echo "$output" | grep "Spatial source:" | grep "W/m" | tail -1 | awk '{print $3}')
	power_bounce=$(echo "$output" | grep "Bounce-averaged:" | grep "W/m" | tail -1 | awk '{print $2}')
	power_midmap=$(echo "$output" | grep "Midplane-mapped:" | grep "W/m" | tail -1 | awk '{print $2}')
	
	# Extract M0 fluxes
	M0_spatial=$(echo "$output" | grep "Spatial source:" | grep "s^-1" | tail -1 | awk '{print $3}')
	M0_contour=$(echo "$output" | grep "Contour integral:" | tail -1 | awk '{print $3}')
	M0_orbit=$(echo "$output" | grep "Orbit-averaged:" | grep "s^-1" | tail -1 | awk '{print $2}')
	
	# Extract average energies - second to last field ($(NF-1)) is the numeric value before "eV"
	E_avg_spatial=$(echo "$output" | grep "Spatial source:" | grep "eV$" | tail -1 | awk '{print $(NF-1)}')
	E_avg_bounce=$(echo "$output" | grep "Bounce-averaged:" | grep "eV$" | tail -1 | awk '{print $(NF-1)}')
	E_avg_midmap=$(echo "$output" | grep "Midplane-mapped:" | grep "eV$" | tail -1 | awk '{print $(NF-1)}')
	
	# Store optimized parameters in arrays
	E_beam_midmap_values+=("$E_beam_opt")
	T_beam_midmap_values+=("$T_beam_opt")
	gamma_midmap_values+=("$gamma_midmap_opt")
	
	echo "Optimized orbit-averaged source parameters:"
	echo "  E_beam_midmap   = $E_beam_opt eV"
	echo "  T_beam_midmap   = $T_beam_opt eV"
	echo "  gamma_midmap    = $gamma_midmap_opt"
	echo ""
	echo "Source comparisons:"
	echo "  Total power deposited (W/m^3):"
	echo "    Spatial source:   $power_spatial"
	echo "    Contour integral:  $power_bounce"
	echo "    Midplane-mapped:  $power_midmap"
	echo ""
	echo "  Total M0 fluxes (s^-1 m^-3):"
	echo "    Spatial source:   $M0_spatial"
	echo "    Contour integral: $M0_contour"
	echo "    Orbit-averaged:   $M0_orbit"
	echo ""
	echo "  Average energy per particle (eV):"
	echo "    Spatial source:   $E_avg_spatial"
	echo "    Contour integral:  $E_avg_bounce"
	echo "    Midplane-mapped:  $E_avg_midmap"
	echo ""
	echo "=========================================="
	echo ""
done

rm temp.py
# Print condensed summary for copy-run.sh
echo ""
echo "=============================================="
echo "COPY-PASTE READY FORMAT"
echo "=============================================="
echo ""
echo "E_beam_midmap_values=(${E_beam_midmap_values[@]})"
echo "T_beam_midmap_values=(${T_beam_midmap_values[@]})"
echo "gamma_midmap_values=(${gamma_midmap_values[@]})"
echo ""
echo "=============================================="
echo "Individual values:"
echo "=============================================="
for i in "${!R_values[@]}"; do
  echo "R=${R_values[$i]}: E_beam=${E_beam_midmap_values[$i]} eV, T_beam=${T_beam_midmap_values[$i]} eV, gamma=${gamma_midmap_values[$i]}"
done
echo "=============================================="
