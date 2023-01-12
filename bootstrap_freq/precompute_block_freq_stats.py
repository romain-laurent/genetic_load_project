#!/usr/bin/env python3

import sys
import gzip

print('Reading score classes...', file=sys.stderr)
f = open('../CADD_scores_computations/scores_with_classes.txt')
score_classes = set()
next(f)                         # header
dico_snp = {}
for line in f :
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    gerp_class, cadd_class = line[-2], line[-1]
    if (chrom, pos) in dico_snp :
        print('Problem multiple classes...', file=sys.stderr)
    score_classes.add(gerp_class)
    score_classes.add(cadd_class)
    dico_snp[(chrom, pos)] = [0, gerp_class, cadd_class]
f.close()

print('Reading block definitions...', file=sys.stderr)
f = open('../bootstrap_blocks_definition.txt')
next(f)
block_ids = set()
for line in f :
    line = line.split()
    block, chrom, pos = int(line[0]), int(line[1]), int(line[2])
    block_ids.add(block)
    dico_snp[(chrom,pos)][0] = block
f.close()

print('Reading genotypes and counting alleles...', file=sys.stderr)
f = gzip.open('../../annotation/refined_86_samples_ref_is_ancestral.vcf.gz')
counts = {}
dico_idxs = {'0/0':0, '1/1':2, '0/1':1, '1/0':1}
for line in f :
    line = line.decode('utf-8')
    if line.startswith('##') :
        continue
    if line.startswith('#CHROM') :
        line = line.split()
        indiv_names = line[9:]
        for name in indiv_names :
            counts[name] = {}
            for block_id in block_ids :
                counts[name][block_id] = {}
                for score_class in score_classes :
                    counts[name][block_id][score_class] = [0,0,0]
        continue
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    coords = (chrom,pos)
    genos = line[9:]
    for idx in range(len(genos)) :
        indiv_name = indiv_names[idx]
        geno = genos[idx]
        idx_to_increase = dico_idxs[geno]
        idx_block = dico_snp[coords][0]
        score_class1 = dico_snp[coords][1]
        score_class2 = dico_snp[coords][2]
        counts[indiv_name][idx_block][score_class1][idx_to_increase] += 1
        counts[indiv_name][idx_block][score_class2][idx_to_increase] += 1
    
f.close()

print('indiv\tblock\tscore.class\thom.anc\thet\thom.der')
for indiv in indiv_names :
    for block in block_ids :
        for score in score_classes :
            new_line = [indiv, str(block), score]
            for i in range(3) :
                new_line.append(str(counts[indiv][block][score][i]))
            print('\t'.join(new_line))
            
