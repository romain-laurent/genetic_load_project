#!/usr/bin/env python3

import sys
import random

pop1 = sys.argv[1]
pop2 = sys.argv[2]
nb_bootstrap = 1000

print('Reading precomputed stats...', file=sys.stderr)
dico_blocks = {}
all_blocks = set()
all_scores = set()
f = open(f'precomputed_Rab_{pop1}_{pop2}.txt')
next(f)                         # header
for line in f :
    line = line.split()
    idx_block, score, s1, s2 = int(line[0]), line[1], float(line[2]), float(line[3])
    all_blocks.add(idx_block)
    all_scores.add(score)
    if idx_block not in dico_blocks :
        dico_blocks[idx_block] = {}
    dico_blocks[idx_block][score] = [s1, s2]
f.close()


all_blocks = sorted(list(all_blocks))
all_scores = sorted(list(all_scores))

def compute_stat(dico_blocks, wanted_blocks, all_scores) :
    stats = {}
    for score in all_scores :
        stats[score] = [0,0]
    for block in wanted_blocks :
        for score in all_scores :
            stats[score][0] += dico_blocks[block][score][0]
            stats[score][1] += dico_blocks[block][score][1]
    for score in all_scores :
        stats[score] = stats[score][0] / stats[score][1]
    return stats

# compute stat on real data
stats = compute_stat(dico_blocks, all_blocks, all_scores)
f = open(f'bootstrap_Rab_{pop1}_{pop2}.txt', 'w')
new_line = ['idx'] + all_scores
print('\t'.join(new_line), file=f)
new_line = ['0']
for score in all_scores :
    new_line.append(str(stats[score]))
print('\t'.join(new_line), file=f)

# perform bootstrap
print('Performing bootstrap...', file=sys.stderr)
for idx_boot in range(nb_bootstrap) :
    new_line = [str(idx_boot+1)]
    wanted_blocks = random.choices(all_blocks, k=len(all_blocks))
    stats = compute_stat(dico_blocks, wanted_blocks, all_scores)
    for score in all_scores :
        new_line.append(str(stats[score]))
    print('\t'.join(new_line), file=f)

f.close()
