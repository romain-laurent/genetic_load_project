#!/usr/bin/env python3

import random

# use that to prevent choosing a random genotype each type
# we have a list of already randomize genotypes to use when needed
nb_random = 100000
random_hets = []

f_in = open('../../annotation/refined_86_samples_ref_is_ancestral.vcf')
f_out = open('randomized.vcf', 'w')

for line in f_in :
    if line.startswith('#') :
        print(line.strip(), file=f_out)
        continue
    line = line.split()
    for i in range(9, len(line)) :
        # if we find an heterozygous genotype, we randomize it
        if line[i] == '0/1' or line[i] == '1/0' :
            if len(random_hets) == 0 : # no more random genotypes -> create some new
                random_hets = random.choices(['1/0','0/1'], k=nb_random)
            # set the new randomized genotype
            line[i] = random_hets.pop()
    print('\t'.join(line), file=f_out)
    

f_in.close()
f_out.close()
