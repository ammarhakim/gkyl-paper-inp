echo "Generating data for confinement time plot..."
echo ""

# Function to extract maximum value from pgkyl info output
extract_max() {
    grep "Maximum:" | awk '{print $3}'
}

# R=3
cd R-scan/R-3
M0_R3=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R3=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R3=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=5
cd R-scan/R-5
M0_R5=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R5=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R5=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=10
cd R-scan/R-10
M0_R10=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R10=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R10=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=15
cd R-scan/R-15
M0_R15=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R15=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R15=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=22
cd R-scan/R-22
M0_R22=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R22=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R22=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

cd R-scan/R-50
M0_R50=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R50=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R50=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
M0mid_R50=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

cd ../stellar-lorentzian1x-orbit-average
M0_R32=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_65.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
M0S_R32=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info 2>&1 | extract_max)
nu_R32=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../stellar-lorentzian1x-orbit-average-R-scan

# Print Python-ready output
echo "# Copy the lines below to Python:"
echo ""
echo "intM0dx_R3 = $M0_R3"
echo "intM0dx_R5 = $M0_R5"
echo "intM0dx_R10 = $M0_R10"
echo "intM0dx_R15 = $M0_R15"
echo "intM0dx_R22 = $M0_R22"
echo "intM0dx_R32 = $M0_R32"
echo "intM0dx_R50 = $M0_R50"
echo ""
echo "intM0Sdx_R3 = $M0S_R3"
echo "intM0Sdx_R5 = $M0S_R5"
echo "intM0Sdx_R10 = $M0S_R10"
echo "intM0Sdx_R15 = $M0S_R15"
echo "intM0Sdx_R22 = $M0S_R22"
echo "intM0Sdx_R32 = $M0S_R32"
echo "intM0Sdx_R50 = $M0S_R50"
echo ""
echo "nu_ii_R3 = $nu_R3"
echo "nu_ii_R5 = $nu_R5"
echo "nu_ii_R10 = $nu_R10"
echo "nu_ii_R15 = $nu_R15"
echo "nu_ii_R22 = $nu_R22"
echo "nu_ii_R32 = $nu_R32"
echo "nu_ii_R50 = $nu_R50"