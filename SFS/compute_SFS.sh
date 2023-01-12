#!/bin/bash

max_indiv=19			# number of samples in TJE pop

pops=(CTA LPO TUR TJE)
for pop in ${pops[@]}
do
    echo $pop
    vcftools --vcf ../../annotation/refined_86_samples_ref_is_ancestral.vcf --keep ../samples_${pop}.txt --max-indv $max_indiv --counts2 --out SFS_${pop}_subsampled &
done
