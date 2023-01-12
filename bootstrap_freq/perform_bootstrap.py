#!/usr/bin/env python3

import sys
import random
nb_bootstrap = 1000

pop1 = sys.argv[1]
pop2 = sys.argv[2]

print('Reading samples IDs...', file=sys.stderr)
dico_pops = {}
all_indivs = set()
for pop in [pop1, pop2] :
    dico_pops[pop] = set()
    f = open(f'../samples_{pop}.txt')
    for line in f :
        dico_pops[pop].add(line.strip())
        all_indivs.add(line.strip())
    f.close()


print('Reading precomputed stats...', file=sys.stderr)
dico_counts = {}
f = open('precomputed_freq_stats.txt')
next(f)
all_blocks = set()
all_scores = set()
for line in f :
    line = line.split()
    indiv = line[0]
    if indiv not in all_indivs : # indiv we are not interested in
        continue
    block = int(line[1])
    score = line[2]
    hom_anc, het, hom_der = int(line[3]), int(line[4]), int(line[5])

    all_blocks.add(block)
    all_scores.add(score)

    if indiv not in dico_counts :
        dico_counts[indiv] = {}
    if block not in dico_counts[indiv] :
        dico_counts[indiv][block] = {}
    if score in dico_counts[indiv][block] :
        print('Problem...', file=sys.stderr)

    dico_counts[indiv][block][score] = [hom_anc,het,hom_der]

f.close()

    
def compute_final_stats(counts, blocks, pops, all_scores) :
    # initialise counts per population
    final_counts = {}
    for pop in pops :
        final_counts[pop] = {}
        for score in all_scores :
            final_counts[pop][score] = [0,0,0,0]

    # add over blocks for each pop
    for pop, indivs in pops.items() :
        for indiv in indivs :
            for block in blocks :
                for score in all_scores :
                    # homozygous ancestral
                    final_counts[pop][score][0] += counts[indiv][block][score][0]
                    # all ancestral
                    final_counts[pop][score][1] += 2*counts[indiv][block][score][0] + counts[indiv][block][score][1]
                    # homozygous derived
                    final_counts[pop][score][2] += counts[indiv][block][score][2]
                    # all derived
                    final_counts[pop][score][3] += 2*counts[indiv][block][score][2] + counts[indiv][block][score][1]

    # divide by sample size (i.e. compute mean over individuals)
    for pop in final_counts :
        nb_indiv = len(pops[pop])
        for score in all_scores :
            for i in range(4) :
                final_counts[pop][score][i] /= nb_indiv
    
    # compute ratios
    ratios = {}
    pop1, pop2 = list(pops.keys())
    for score in all_scores :
        ratios[score] = [-1, -1, -1, -1]
        for i in range(4) :
            ratios[score][i] = final_counts[pop1][score][i] / final_counts[pop2][score][i]
    return ratios

def print_ratios(ratios, all_scores, idx_boot) :
    for score in all_scores :
        new_line = [str(idx_boot), score]
        for i in range(4) :
            new_line.append(str(ratios[score][i]))
        print('\t'.join(new_line))
    
print('idx\tscore.class\thom.anc\tall.anc\thom.der\tall.der')

# compute real stats
all_blocks = list(all_blocks)
ratios = compute_final_stats(dico_counts, all_blocks, dico_pops, all_scores)
print_ratios(ratios, all_scores, 0)

# compute bootstrap
print('Computing bootstraps...', file=sys.stderr)
for idx_boot in range(nb_bootstrap) :
    wanted_blocks = random.choices(all_blocks, k=len(all_blocks))
    ratios = compute_final_stats(dico_counts, wanted_blocks, dico_pops, all_scores)
    print_ratios(ratios, all_scores, idx_boot+1)

