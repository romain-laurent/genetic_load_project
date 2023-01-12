#!/usr/bin/env python3

import sys

# get indivs ROH localisations
print('Reading ROH positions for each indiv...', file=sys.stderr)
f = open('../garlic/ROH_classes_2_groups.txt')
next(f)
dico_indivs = {}
for line in f :
    line = line.split()
    pop = line[0]
    indiv = line[1]
    chrom, start, stop = int(line[2]), int(line[3]), int(line[4])
    roh_class = line[-1]
    if pop not in dico_indivs :
        dico_indivs[pop] = {}
    if indiv not in dico_indivs[pop] :
        dico_indivs[pop][indiv] = {}
    if chrom not in dico_indivs[pop][indiv] :
        dico_indivs[pop][indiv][chrom] = []
    dico_indivs[pop][indiv][chrom].append((start, stop, roh_class))
f.close()

# sanity check that ROHs are sorted ('find_roh_type' relies on that assumption)
print('Checking that ROHs are sorted...', file=sys.stderr)
for pop in dico_indivs :
    for indiv in dico_indivs[pop] :
        for chrom in dico_indivs[pop][indiv] :
            prev_pos = -1
            for roh in dico_indivs[pop][indiv][chrom] :
                if roh[0] <= prev_pos :
                    print('Problem: unsorted ROHs...', file=sys.stderr)
                prev_pos = roh[1]

print('Reading SNPs annotations...', file=sys.stderr)
f = open('../CADD_scores_computations/scores_with_classes.txt')
next(f)
annots = {}
all_annots = set()
for line in f :
    line = line.split()
    chrom, pos, annot = int(line[0]), int(line[1]), line[-1]
    if (chrom, pos) in annots :
        print('Problem', file=sys.stderr)
    all_annots.add(annot)
    annots[(chrom, pos)] = annot
f.close()

def find_roh_type(pos, rohs) :
    for roh in rohs :
        if pos > roh[1] :       # pos of interest after end of this ROH
            continue
        if pos < roh[0] :       # pos of interest before this ROH -> we passed it
            break
        if roh[0] <= pos and pos <= roh[1] : # if pos in ROH
            return roh[-1]
    return 'N'

print('Reading genotypes and computing stuff...', file=sys.stderr)
f = open('../../annotation/refined_86_samples_ref_is_ancestral.vcf')
dico_results = {}
for line in f :
    if line.startswith('##') :
        continue
    if line.startswith('#CHROM') : # header -> get indivs idxs
        line = line.split()
        indiv_idxs = {}
        for pop in dico_indivs :
            for indiv in dico_indivs[pop] :
                indiv_idxs[indiv] = line.index(indiv)
                dico_results[indiv] = {}
                for roh_type in ('N','A','B','C') :
                    dico_results[indiv][roh_type] =  {}
                    for annot in all_annots :
                        dico_results[indiv][roh_type][annot] = [0,0,0]
        continue
    # here we have genotypes
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    annot = annots[(chrom,pos)] # get SNP annotation
    # for each indiv
    for pop in dico_indivs :
        for indiv in dico_indivs[pop] :
            genotype = line[indiv_idxs[indiv]] # get genotype
            # find ROH type
            roh_type = find_roh_type(pos, dico_indivs[pop][indiv][chrom])
            # increase relevant variables
            dico_results[indiv][roh_type][annot][0] += 1
            if genotype == '1/1' :
                dico_results[indiv][roh_type][annot][1] += 1
            if genotype == '0/1' or genotype == '1/0' :
                dico_results[indiv][roh_type][annot][2] += 1
        
f.close()
new_line = ['pop','indiv','roh.class','mut.class','hom','het','tot']
print('\t'.join(new_line))
# print results
for pop in dico_indivs :
    for indiv in dico_indivs[pop] :
        for roh_class in ('N','A','B','C') :
            for annot in all_annots :
                new_line = [pop, indiv, roh_class, annot]
                new_line.append(str(dico_results[indiv][roh_class][annot][1]))
                new_line.append(str(dico_results[indiv][roh_class][annot][2]))
                new_line.append(str(dico_results[indiv][roh_class][annot][0]))
                print('\t'.join(new_line))
