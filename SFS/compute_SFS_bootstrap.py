#!/usr/bin/env python3

import random

nb_sample_wanted = 19
nb_bootstrap = 100

pops = ['CTA','LPO','TUR','TJE']

# get sample ids for each pop
dico_samples = {}
for pop in pops :
    f = open(f'../samples_{pop}.txt')
    lines = f.readlines()
    f.close()
    dico_samples[pop] = [line.strip() for line in lines]


# generate datast to keep results
nb_possible_derived_alleles = 2*nb_sample_wanted - 1
dico_results = {}
for pop in pops :
    dico_results[pop] = []
    for i in range(nb_bootstrap) :
        tmp = [0] * nb_possible_derived_alleles
        dico_results[pop].append(tmp)

# generate subsamples of indiviuals for each bootstrap
wanted_indivs = {}
for pop in pops :
    wanted_indivs[pop] = []
    for i in range(nb_bootstrap) :
        tmp = random.choices(dico_samples[pop], k=nb_sample_wanted)
        wanted_indivs[pop].append(tmp)

# here we go
f = open('../../annotation/refined_86_samples_ref_is_ancestral.vcf')
idx_snp = 0
for line in f :
    idx_snp += 1
    if (idx_snp % 1000 == 0) :
        print(idx_snp)
    if line.startswith('##') :
        continue
    # header -> find samples columns index in VCF file
    if line.startswith('#') :
        dico_idxs = {}
        line = line.split()
        for pop in pops :
            for sample in dico_samples[pop] :
                dico_idxs[sample] = line.index(sample)
        continue
    line = line.split()
    # for each pop
    for pop in pops :
        # for each bootstrap
        for idx_boot in range(nb_bootstrap) :
            tmp_count = 0
            # we count the number of derived alleles for this SNP in this bootstrap
            for sample in wanted_indivs[pop][idx_boot] :
                tmp_count += line[dico_idxs[sample]].count('1')
            # store the result if not monomorphic
            if tmp_count > 0 and tmp_count < 2*nb_sample_wanted :
                dico_results[pop][idx_boot][tmp_count-1] += 1
f.close()

for pop in pops :
    f = open(f'result_bootstrap_{pop}.txt','w')
    for idx_boot in range(nb_bootstrap) :
        new_line = [str(i) for i in dico_results[pop][idx_boot]]
        print('\t'.join(new_line), file=f)
    f.close()
