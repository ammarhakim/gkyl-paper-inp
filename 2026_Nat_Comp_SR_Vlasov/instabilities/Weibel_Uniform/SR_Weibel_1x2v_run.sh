#!/bin/bash
##### Scan the parameter space for the instability growth rate
##### run with: bash SCAN_GROWTH_RATE.sh

##### Make a file to save all the output text
rm -f BACKUP_SR_WEIBEL_1X2V/GROWTH_RATE.txt
rm -f BACKUP_SR_WEIBEL_1X2V/out.txt
rm -f BACKUP_SR_WEIBEL_1X2V/Growth_Plot_*
touch "BACKUP_SR_WEIBEL_1X2V/GROWTH_RATE.txt"
touch "BACKUP_SR_WEIBEL_1X2V/out.txt"

##### Activate postgkyl
conda activate pgkyl

##### Integer counter for file names: 1, 2, 3, ...
i=1

##### Select a range of k values (doubles)
for k in $(seq 0.1 0.5 2.6)
do
    ##### Format k nicely (e.g. 0.0, 0.1, 0.2, ...)
    replace_str=$(printf "%.1f" "$k")

    ##### Optionally, padded integer index if you like 001, 002, ...
    idx=$i           # or: idx=$(printf "%03d" "$i")

    scan_var="kx"
    k_init="0.4"

    echo "Loop $idx: replacing $scan_var = $k_init with $scan_var = $replace_str using sed"

    sed "s/$scan_var = $k_init/$scan_var = $replace_str/" \
        TEST_sr_weibel_1x2v.c > vlasov/creg/rt_scan_weibel_sr_1x2v.c

    ##### Compile and Run the code
    make vlasov-regression -j 32
    ./build/vlasov/creg/rt_scan_weibel_sr_1x2v #####-g

    ##### When done, run pgkyl, save the image for the growth rate
    pgkyl vlasov_sr_weibel_1x2v-field-energy.gkyl sel -c5 growth -d plot -f0 --logy --no-show \
        --saveas "BACKUP_SR_WEIBEL_1X2V/Growth_Plot_${scan_var}_${idx}.png" \
        >> BACKUP_SR_WEIBEL_1X2V/out.txt

    echo "$scan_var = $replace_str (index $idx)" >> BACKUP_SR_WEIBEL_1X2V/GROWTH_RATE.txt
    tail -n 1 BACKUP_SR_WEIBEL_1X2V/out.txt >> BACKUP_SR_WEIBEL_1X2V/GROWTH_RATE.txt

    ##### Create a copy of the file so it doesn't need to be run again:
    cp vlasov_sr_weibel_1x2v-field-energy.gkyl \
       "BACKUP_SR_WEIBEL_1X2V/vlasov_sr_weibel_1x2v-field-energy_${idx}.gkyl"

    ##### Increment integer index
    i=$((i + 1))

##### End of the Loop/File
done
