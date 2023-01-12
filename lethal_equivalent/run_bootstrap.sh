#!/bin/bash

neutral=VEP_SS
non_neutral=VEP_LOF

#./compute_bootstrap.py $neutral $non_neutral > bootstrap_${neutral}_${non_neutral}.txt &


neutral=CADD_inf_10
non_neutral=CADD_sup_28

./compute_bootstrap.py $neutral $non_neutral > bootstrap_${neutral}_${non_neutral}.txt  2> bootstrap_${neutral}_${non_neutral}.err &

