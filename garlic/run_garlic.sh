#!/bin/bash

pops=(CTA LPO TUR TJE)

for pop in ${pops[@]}
do
    ./garlic \
        --weighted \
        --auto-overlap-frac \
        --auto-winsize \
        --error 0.001 \
        --tfam ready_garlic_$pop.tfam \
        --tped ready_garlic_$pop.tped \
        --cm \
        --build hg19 \
        --map genetic_map_$pop.txt \
        --out results_weighted_scheme/ROH_$pop \
        --threads 30


done

