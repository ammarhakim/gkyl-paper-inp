#!/bin/bash

frNum=0
elem1=Li0
elem2=Li1
elem3=Li2
elem4=Ar8
elem5=Ar9
elem6=Ar11
echo "#    Te    Li0_Max  Li0_NM  Li1_Max  Li1_NM  Li2_Max  Li2_NM  Ar8_Max  Ar8_NM  Ar9_Max  Ar9_NM  Ar11_Max  Ar11_NM" >all_emissivities.txt
#for i in {0..198};
#do
#    ./non_maxwellian_compare -s 1 -w $i
#    j=$((i+1))
#    mv *.json Te${j}
#    mv *.gkyl Te${j}
#done

#conda activate pgkyl
for i in {1..199};
do
    simName=Te${i}/non_maxwellian_elc
    echo $i "  " $(pgkyl "$simName-elc_M0_$frNum.gkyl" -t n "$simName-elc_M1_$frNum.gkyl" -t m1 "$simName-elc_M2_$frNum.gkyl" -t m2 interp -b ms -p1 ev 'm2 m1 m1 * n / - n / 9.1672621e-31 * 3 /  1.6021766e-19 /' info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem1}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem1}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem2}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem2}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem3}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem3}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem4}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem4}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem5}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem5}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_maxwellian_${elem6}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1) $(pgkyl "$simName-elc_radiation_emissivity_${elem6}_0.gkyl" interp -b ms -p 1 info | tail -n 6 | awk '{print $3}' | head -n 1)>> all_emissivities.txt
done
