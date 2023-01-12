#!/usr/bin/env python3

import sys
import gzip

pop1 = sys.argv[1]
pop2 = sys.argv[2]

# read samples
pops = (pop1, pop2)
dico_pops = {}
print('Reading sample infos...', file=sys.stderr)
for pop in pops :
    dico_pops[pop] = set()
    f = open(f'../samples_{pop}.txt')
    for line in f :
        dico_pops[pop].add(line.strip())
    f.close()

# read blocks
print('Reading blocks infos...', file=sys.stderr)
dico_snps = {}
f = open('../bootstrap_blocks_definition.txt')
next(f)                         # header
all_blocks = set()
for line in f :
    line = line.split()
    idx, chrom, pos = int(line[0]), int(line[1]), int(line[2])
    dico_snps[(chrom, pos)] = idx
    all_blocks.add(idx)
f.close()

print('Reading score classes...', file=sys.stderr)
f = open('../CADD_scores_computations/scores_with_classes.txt')
next(f)
all_scores = set()
for line in f :
    line = line.split()
    coords = (int(line[0]), int(line[1]))
    block = dico_snps[coords]
    gerp, cadd = line[-2], line[-1]
    dico_snps[coords] = (block, gerp, cadd)
    all_scores.add(gerp)
    all_scores.add(cadd)
f.close()

def compute_freqs(genos, dico_pops, dico_idxs) :
    freqs = {}
    for pop in dico_pops :
        freq = 0
        for indiv in dico_pops[pop] :
            geno = genos[dico_idxs[indiv]]
            if geno == '1/1' :
                freq += 2
            elif geno == '1/0' or geno == '0/1' :
                freq += 1
        freqs[pop] = freq / (2*len(dico_pops[pop]))
    return freqs
            
print('Reading genotypes and computing stats...', file=sys.stderr)
all_stats = {}
for block in all_blocks :
    all_stats[block] = {}
    for score in all_scores :
        all_stats[block][score] = [0,0]
        
f = open('../../annotation/refined_86_samples_ref_is_ancestral.vcf')
for line in f :

    if line.startswith('##') :
        continue
    if line.startswith('#CHROM') : # header -> find indivs idxs
        line = line.split()
        dico_idxs = {}
        for pop in pops :
            for indiv in dico_pops[pop] :
                dico_idxs[indiv] = line.index(indiv)
        continue
    line = line.split()
    freqs = compute_freqs(line, dico_pops, dico_idxs)
    stat_1_not_2 = freqs[pop1] * (1-freqs[pop2])
    stat_2_not_1 = freqs[pop2] * (1-freqs[pop1])
    # store stat
    coords = (int(line[0]), int(line[1]))
    block, score1, score2 = dico_snps[coords]
    all_stats[block][score1][0] += stat_1_not_2
    all_stats[block][score1][1] += stat_2_not_1
    all_stats[block][score2][0] += stat_1_not_2
    all_stats[block][score2][1] += stat_2_not_1
f.close()

all_blocks = sorted(list(all_blocks))
all_scores = sorted(list(all_scores))
f = open(f'precomputed_Rab_{pop1}_{pop2}.txt', 'w')
new_line = ['block', 'score', f'L_{pop1}_not_{pop2}', f'L_{pop2}_not_{pop1}']
print('\t'.join(new_line), file=f)
for block in all_blocks :
    for score in all_scores :
        new_line = [str(block), score]
        new_line.append(str(all_stats[block][score][0]))
        new_line.append(str(all_stats[block][score][1]))
        print('\t'.join(new_line), file=f)
f.close()

