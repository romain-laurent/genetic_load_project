#!/bin/bash

pops=(CTA LPO TUR TJE)
for pop in ${pops[@]}
do
    tail --lines=+2 ../LOF/${pop}_Table_analyse_variant_count_VEP_SNP_derived.txt | awk '{print $1"\t"$2}' > variants_LOF_$pop.txt
    vcftools --vcf ../VCFs/${pop}_final.vcf --positions variants_LOF_$pop.txt --recode --out ${pop}_LOF_annotated
done
rm variants_LOF*
rm *log
