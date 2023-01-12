#!/usr/bin/env python3

import sys
import gzip

pop = sys.argv[1]
print(f'Converting to plink format for pop {pop}...', file=sys.stderr)

f = open(f'freqs_{pop}.frq')
alleles = []
next(f)
for line in f :
    line = line.split()
    chrom, a1, a2, name = int(line[0]), line[2], line[3], line[1]
    nb_chrom = int(line[-1])
    alleles.append( (chrom, a1, a2, name) )
f.close()

f_in = gzip.open(f'{pop}.freq.gz')
next(f_in)
f_out = open(f'final_freqs_{pop}.frq', 'w')
print('CHR SNP A1 A2 MAF NCHROBS', file=f_out)
idx_line = 0
for line in f_in :
    line = line.decode('utf8')
    line = line.split()
    chrom = int(line[0][3:])
    name = line[1]
    pos = line[2]
    allele_garlic = line[3]
    freq = float(line[-1])
    if alleles[idx_line][-1] != name :
        print(f'problem name {name}', file=sys.stderr)
        print(f'{idx_line} {alleles[idx_line]}', file=sys.stderr)
        sys.exit(1)
    if alleles[idx_line][0] != chrom :
        print(f'problem chrom {chrom} {pop}', file=sys.stderr)
    a1, a2 = alleles[idx_line][1], alleles[idx_line][2]
    if a1 != allele_garlic and a2 != allele_garlic :
        print(f'problem allele {chrom} {pop} {line[1]}', file=sys.stderr)
        sys.exit()
    if freq > 0.5 :
        if allele_garlic == a1 :
            a1, a2 = a2, a1
        freq = 1 - freq
    new_line = f'{chrom} {name} {a1} {a2} {freq} {nb_chrom}'
    print(new_line, file=f_out)
    idx_line += 1
f_in.close()
f_out.close()
