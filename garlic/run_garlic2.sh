#!/bin/bash

# pops where weighted scheme did not really work
pops=(CTA LPO TUR TJE)

for pop in ${pops[@]}
do
    
    for win_size in {30..100..10}
    do

    ./garlic \
        --weighted \
        --auto-overlap-frac \
        --winsize $win_size \
        --error 0.001 \
        --tfam ready_garlic_$pop.tfam \
        --tped ready_garlic_$pop.tped \
        --cm \
        --build hg19 \
        --map genetic_map_$pop.txt \
        --out results_force_winsize/ROH_${pop}_winsize_$win_size \
        --threads 30
    done

    
done

