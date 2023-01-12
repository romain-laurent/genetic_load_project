#!/usr/bin/env python3

import sys

pops = ['CTA','LPO','TUR','TJE']
dico_indivs = {}
for pop in pops :
    f = open(f'../samples_{pop}.txt')
    dico_indivs[pop] = []
    for line in f :
        line = line.strip()
        dico_indivs[pop].append(line)
    f.close()

annotations = {}
f = open('../VEP/VEP_annotations.txt')
next(f)                         # header
for line in f :
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    annotations[(chrom, pos)] = line[2]
f.close()

def analyse_genos(line, idxs) :
    nb_hom = 0
    nb_alt = 0
    for idx in idxs :
        geno = line[idx]
        if geno == '1/1' :
            nb_hom += 1
            nb_alt += 2
        elif geno == '0/1' or geno == '1/0' :
            nb_alt += 1
    return nb_hom, nb_alt

dico_idxs = {}
dico_snps = {}
f = open('86_samples_ref_ancestral_LOF.recode.vcf')
for line in f :
    if line.startswith('##') :
        continue
    if line.startswith('#') :
        line = line.split()
        for pop in pops :
            dico_idxs[pop] = []
            for indiv in dico_indivs[pop] :
                dico_idxs[pop].append(line.index(indiv))
        continue
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    if (chrom, pos) not in annotations :
        print('Unknown annotations...', file=sys.stderr)
    dico_snps[(chrom,pos)] = {}
    for pop in pops :
        hom_alt, nb_alt = analyse_genos(line, dico_idxs[pop])
        dico_snps[(chrom, pos)][pop] = (hom_alt, nb_alt)
        
f.close()

print('\t'.join(['chrom','pos','annot','pop','hom','alt']))
for coords in dico_snps :
    for pop in pops :
        new_line = [str(coords[0]), str(coords[1])]
        new_line.append(annotations[coords])
        new_line.append(pop)
        new_line.append(str(dico_snps[coords][pop][0]))
        new_line.append(str(dico_snps[coords][pop][1]))
        print('\t'.join(new_line))
