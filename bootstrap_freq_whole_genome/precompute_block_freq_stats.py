#!/usr/bin/env python3

import sys

# get blocks definitions
positions = {}
f = open('../bootstrap_blocks_definition.txt')
next(f)                         # header
for line in f :
    line = line.split()
    idx_block, chrom, pos = [int(i) for i in line]
    if (chrom, pos) in positions :
        print('duplicated position')
    positions[(chrom,pos)] = idx_block
f.close()

f = open('../../annotation/refined_86_samples_ref_is_ancestral.vcf')
counts = {}
for line in f :
    if line.startswith('##') :
        continue
    if line.startswith('#CHROM') :
        indiv_names = line.split()[9:]
        for indiv in indiv_names :
            counts[indiv] = {}
        continue
    line = line.split()
    
    chrom, pos = int(line[0]), int(line[1])
    idx_block = positions[(chrom, pos)]
    genos = line[9:]
    
    for idx_geno in range(len(genos)) :
        indiv = indiv_names[idx_geno]
        if idx_block not in counts[indiv] :
            counts[indiv][idx_block] = [0, 0, 0]
        geno = genos[idx_geno]
        if geno == '0/0' :
            counts[indiv][idx_block][0] += 1
        elif geno == '1/1' :
            counts[indiv][idx_block][1] += 1
        elif geno == '0/1' or geno == '1/0' :
            counts[indiv][idx_block][2] += 1
        else :
            print('Weird genotype',file=sys.stderr)

f.close()

# we counted everything, write the result
print('indiv\tblock\thom.anc\thom.der\thet')
for indiv in indiv_names :
    for i in range(len(counts[indiv])) :
        h_anc, h_der, het = counts[indiv][i+1]
        print(f'{indiv}\t{i+1}\t{h_anc}\t{h_der}\t{het}')
