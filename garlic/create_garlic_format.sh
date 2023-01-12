#!/bin/bash

# WE WORK ON THE FILE WITH AS MUCH VARIANTS AS POSSIBLE, I.E. BEFORE ASSIGNING ANCESTRAL/DERIVED INFO

# need to name the variants to identify problems later...
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ../../annotation/refined_86_samples_no_annots.vcf.gz -Oz -o refined_86_samples_no_annots_with_SNP_ids.vcf.gz

# create one file per population
pops=(CTA LPO TUR TJE)
for pop in ${pops[@]}
do
    echo $pop
#    bcftools view -S ../samples_${pop}.txt refined_86_samples_no_annots_with_SNP_ids.vcf.gz -Ov -o ${pop}_from_86_samples.vcf
    
#    plink --vcf ${pop}_from_86_samples.vcf --cm-map ../Festim/genetic_maps/map_chr@.gmap --recode transpose --out ready_garlic_$pop
    awk '{print $1,$2,$3,$4}' ready_garlic_$pop.tped > genetic_map_$pop.txt
    awk -v pop=$pop '{$1=pop; print}' ready_garlic_$pop.tfam > tmp
    mv tmp ready_garlic_$pop.tfam
done



rm *log *nosex *vcf
