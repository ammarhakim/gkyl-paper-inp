#!/bin/bash
##### Scan the parameter space for the instability growth rate
##### run with: bash SCAN_GROWTH_RATE.sh

##### Make a file to save all the output text
rm -f BACKUP_NU_SR_TWOSTREAM_1X1V/GROWTH_RATE.txt
rm -f BACKUP_NU_SR_TWOSTREAM_1X1V/out.txt
rm -f BACKUP_NU_SR_TWOSTREAM_1X1V/Growth_Plot_*
touch "BACKUP_NU_SR_TWOSTREAM_1X1V/GROWTH_RATE.txt"
touch "BACKUP_NU_SR_TWOSTREAM_1X1V/out.txt"

##### Activate postgkyl
conda activate pgkyl

##### Integer counter for file names: 1, 2, 3, ...
i=1

##### Select a range of k values (doubles)
for k in $(seq 0.03 0.02 0.13)
do
    ##### Format k nicely (e.g. 0.00, 0.10, 0.20, ...)
    replace_str=$(printf "%.2f" "$k")

    ##### Optionally, padded integer index if you like 001, 002, ...
    idx=$i           # or: idx=$(printf "%03d" "$i")

    scan_var="kx"
    k_init="0.02"

    echo "Loop $idx: replacing $scan_var = $k_init with $scan_var = $replace_str using sed"

    sed "s/$scan_var = $k_init/$scan_var = $replace_str/" \
        vlasov/luareg/rt_test_vlasov_sr_nonuniformv_twostream_1x1v.lua > vlasov/luareg/rt_scan_vlasov_sr_nonuniformv_twostream_1x1v.lua

    ##### Compile and Run the code
    make vlasov-regression -j 32
    ./build/gkeyll/gkeyll vlasov/luareg/rt_scan_vlasov_sr_nonuniformv_twostream_1x1v.lua

    ##### When done, run pgkyl, save the image for the growth rate
    pgkyl rt_scan_vlasov_sr_nonuniformv_twostream_1x1v-field-energy.gkyl sel -c0 growth -d plot -f0 --logy --no-show \
        --saveas "BACKUP_NU_SR_TWOSTREAM_1X1V/Growth_Plot_${scan_var}_${idx}.png" \
        >> BACKUP_NU_SR_TWOSTREAM_1X1V/out.txt

    echo "$scan_var = $replace_str (index $idx)" >> BACKUP_NU_SR_TWOSTREAM_1X1V/GROWTH_RATE.txt
    tail -n 1 BACKUP_NU_SR_TWOSTREAM_1X1V/out.txt >> BACKUP_NU_SR_TWOSTREAM_1X1V/GROWTH_RATE.txt

    ##### Create a copy of the file so it doesn't need to be run again:
    cp rt_scan_vlasov_sr_nonuniformv_twostream_1x1v-field-energy.gkyl \
       "BACKUP_NU_SR_TWOSTREAM_1X1V/rt_scan_vlasov_sr_nonuniformv_twostream_1x1v-field-energy_${idx}.gkyl"

    ##### Increment integer index
    i=$((i + 1))

##### End of the Loop/File
done
