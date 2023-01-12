#!/bin/bash

pops=(CTA LPO TUR TJE)

for pop in ${pops[@]}
do
    echo $pop
    # put in tped/tfam format for garlic
    plink --bfile ready_festim_$pop --recode transpose --out for_garlic_$pop
    # use garlic to compute allele freqs with resampling
     ./garlic --tped for_garlic_$pop.tped --tfam for_garlic_$pop.tfam --build hg19 --winsize 2 --error 0.000001 --freq-only --resample 100 --out $pop
    
    # need to use plink as well, to later convert to plink format
     plink --bfile ready_festim_$pop --freq --out freqs_$pop
     ./convert_to_plink_frq.py $pop

done

rm *log *error *nosex 
rm for_garlic*
rm freqs* *gz
