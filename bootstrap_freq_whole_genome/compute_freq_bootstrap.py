#!/usr/bin/env python3

import sys
import random

if len(sys.argv) != 3 :
    print(f'usage: {sys.argv[0]} pop1 pop2', file=sys.stderr)
    sys.exit(1)

pop1 = sys.argv[1]
pop2 = sys.argv[2]

    
# get indiv lists
indivs = {}
for pop in [pop1, pop2] :
    indivs[pop] = set()
    f = open(f'../samples_{pop}.txt')
    for line in f :
        indivs[pop].add(line.strip())
    f.close()

all_indivs = indivs[pop1].union(indivs[pop2])

# read block infos for indivs we are interested in
blocks = {}
all_blocks = set()
f = open('counts_per_block.txt')
next(f)                        # header
for line in f :
    line = line.split()
    indiv = line[0]
    if indiv not in all_indivs :
        continue
    if indiv not in blocks :
        blocks[indiv] = {}
    block_id, hom_anc, hom_der, het = int(line[1]), int(line[2]), int(line[3]), int(line[4])
    all_blocks.add(block_id)
    if block_id in blocks[indiv] :
        print('big problem', file=sys.stderr)
    blocks[indiv][block_id] = (hom_anc, hom_der, het)

f.close()

# compute the ratio (per pop) of mean number of derived allele per indiv
# we actually compute two ratios: one with only derived homozygotes, the other with derived hom + het
def compute_stat(indivs, blocks, wanted_blocks, pop1, pop2) :
    nb_derived = {}
    for indiv in blocks :
        nb_derived[indiv] = [0, 0, 0, 0]
        for block_id in wanted_blocks :
            n_anc_hom = blocks[indiv][block_id][0] * 2
            n_anc_tot = n_anc_hom + blocks[indiv][block_id][2]
            n_deriv_hom = blocks[indiv][block_id][1] * 2
            n_deriv_tot = n_deriv_hom + blocks[indiv][block_id][2]
            nb_derived[indiv][0] += n_anc_hom
            nb_derived[indiv][1] += n_anc_tot
            nb_derived[indiv][2] += n_deriv_hom
            nb_derived[indiv][3] += n_deriv_tot
    means = {pop1:[0,0,0,0], pop2:[0,0,0,0]}
    for pop in [pop1, pop2] :
        for indiv in indivs[pop] :
            means[pop][0] += nb_derived[indiv][0]
            means[pop][1] += nb_derived[indiv][1]
            means[pop][2] += nb_derived[indiv][2]
            means[pop][3] += nb_derived[indiv][3]
        means[pop][0] /= len(indivs[pop])
        means[pop][1] /= len(indivs[pop])
        means[pop][2] /= len(indivs[pop])
        means[pop][3] /= len(indivs[pop])
    ratios = [0, 0, 0, 0]
    ratios[0] = means[pop1][0] / means[pop2][0]
    ratios[1] = means[pop1][1] / means[pop2][1]
    ratios[2] = means[pop1][2] / means[pop2][2]
    ratios[3] = means[pop1][3] / means[pop2][3]
    return ratios
            

# compute stat on the whole genome
all_blocks = list(all_blocks)
real_stat = compute_stat(indivs, blocks, all_blocks, pop1, pop2)

print('idx\thom_anc\tall_anc\thom_der\tall_der')
print(f'0\t{real_stat[0]}\t{real_stat[1]}\t{real_stat[2]}\t{real_stat[3]}')
# compute bootstraps
nb_bootstrap = 1000
for idx_boot in range(nb_bootstrap) :
    wanted_idxs = random.choices(all_blocks, k=len(all_blocks))
    tmp_stat = compute_stat(indivs, blocks, wanted_idxs, pop1, pop2)
    print(f'{idx_boot+1}\t{tmp_stat[0]}\t{tmp_stat[1]}\t{tmp_stat[2]}\t{tmp_stat[3]}')
