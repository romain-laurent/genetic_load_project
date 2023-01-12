#!/bin/bash

pop=$1

../../FSuite/FSuite_1.0.4/FSuite/fsuite.pl --create-submaps --map ready_festim_$pop.bim --n-submaps 100 --hotspots hg19 --sub-folder submaps-$pop
../../FSuite/FSuite_1.0.4/FSuite/fsuite.pl --FEstim --bfile ready_festim_$pop --freq final_freqs_$pop.frq --sub-folder submaps-$pop --mating-type --out festim-$pop

awk '{print $1, $2}' festim-$pop.summary | tail --lines=+2 > festim-$pop.list
../../FSuite/FSuite_1.0.4/FSuite/fsuite.pl --FLOD --bfile ready_festim_$pop --freq final_freqs_$pop.frq --sub-folder submaps-$pop --list-id festim-$pop --FLOD-folder flod-$pop
