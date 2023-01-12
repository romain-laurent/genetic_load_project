#!/bin/bash

# WE WORK ON THE FILE WITH AS MUCH VARIANTS AS POSSIBLE, I.E. BEFORE ASSIGNING ANCESTRAL/DERIVED INFO

# need to name the variants to identify problems later...
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ../../annotation/refined_86_samples_no_annots.vcf.gz -Oz -o refined_86_samples_no_annots_with_SNP_ids.vcf.gz

pops=(CTA LPO TUR TJE)
for pop in ${pops[@]}
do
    echo $pop
    # Festim needs to be run pop by pop -> create one file for each pop
    # also remove markers that become totally monomorphic because they are totally useless for Festim anyway
#    bcftools view -S ../samples_${pop}.txt --min-ac 1:minor refined_86_samples_no_annots_with_SNP_ids.vcf.gz -Ov -o ${pop}_no_monomorphs_from_86_samples.vcf
    

    
    plink --vcf ${pop}_no_monomorphs_from_86_samples.vcf --cm-map genetic_maps/map_chr@.gmap --make-bed --out ready_festim_$pop
    rm *log *nosex
    rm ${pop}_no_monomorphs_from_86_samples.vcf
done




