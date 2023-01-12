#!/usr/bin/env python3

import sys

pops = ['CTA','LPO','TUR','TJE']


print('Reading samples infos...', file=sys.stderr)
dico_indivs = {}
# read samples
for pop in pops :
    dico_indivs[pop] = []
    f = open(f'../samples_{pop}.txt')
    for line in f :
        dico_indivs[pop].append(line.strip())
    f.close()

print('Reading annotations infos...', file=sys.stderr)
dico_annots = {}
f = open('../VEP/VEP_annotations.txt')
next(f)                         # header
for line in f :
    line = line.split()
    coords = (int(line[0]), int(line[1]))
    if coords in dico_annots :
        print('Problem: duplicated annotation...', file=sys.stderr)
    dico_annots[coords] = line[2]
f.close()

print('Reading ROH infos...', file=sys.stderr)
f = open('../garlic/ROH_classes_2_groups.txt')
next(f)                         # header
dico_roh = {}
all_roh_classes = set()
for line in f :
    line = line.split()
    indiv = line[1]
    chrom, start, stop = int(line[2]), int(line[3]), int(line[4])
    roh_class = line[-1]
    all_roh_classes.add(roh_class)
    roh = (start, stop, roh_class)
    if indiv not in dico_roh :
        dico_roh[indiv] = {}
    if chrom not in dico_roh[indiv] :
        dico_roh[indiv][chrom] = []
    dico_roh[indiv][chrom].append(roh)
f.close()

def find_roh_type(pos, rohs) :
    for roh in rohs :
        if roh[0] > pos :       # assumes ROH are sorted in ascending order
            return 'N'
        if pos >= roh[0] and pos <= roh[1] :
            return roh[2]
    return 'N'

print('Counting variants...', file=sys.stderr)
all_types = set()
dico_counts = {}
for coords in dico_annots :
    chrom, pos = coords
    for indiv in dico_roh :
        roh_type = find_roh_type(pos, dico_roh[indiv][chrom])
        annot = dico_annots[coords]
        to_add = f'{annot}_{roh_type}'
        all_types.add(to_add)
        if indiv not in dico_counts :
            dico_counts[indiv] = {}
        if to_add not in dico_counts[indiv] :
            dico_counts[indiv][to_add] = 0
        dico_counts[indiv][to_add] += 1

# print result
new_line = ['pop','indiv']
for tmp in all_types :
    new_line.append(tmp)
print('\t'.join(new_line))
for pop in pops :
    for indiv in dico_indivs[pop] :
        new_line = [pop, indiv]
        for tmp in all_types :
            if tmp not in dico_counts[indiv] :
                new_line.append('0')
            else :
                new_line.append(str(dico_counts[indiv][tmp]))
        print('\t'.join(new_line))
