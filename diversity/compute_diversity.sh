#!/bin/bash

pops=( CTA LPO TUR TJE )

winsize=1000

for pop in ${pops[@]}
do
    #vcftools --vcf ../../annotation/refined_86_samples_ref_is_ancestral.vcf --keep ../samples_$pop.txt --window-pi ${winsize}000 --out diversity_${pop}_win${winsize}k &
    vcftools --vcf randomized.vcf --keep ../samples_$pop.txt --window-pi ${winsize}000 --out diversity_randomized_hets_${pop}_win${winsize}k &
done
