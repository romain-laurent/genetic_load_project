#!/usr/bin/env python3

f = open('../garlic/ROH_classes_2_groups.txt')
dico = {}
all_roh_classes = set()
next(f)
for line in f :
    line = line.split()
    pop, indiv = line[0], line[1]
    roh_class = line[-1]
    if roh_class not in all_roh_classes :
        all_roh_classes.add(roh_class)
    length = int(line[4]) - int(line[3])
    if pop not in dico :
        dico[pop] = {}
    if indiv not in dico[pop] :
        dico[pop][indiv] = {}
    if roh_class not in dico[pop][indiv] :
        dico[pop][indiv][roh_class] = 0
    dico[pop][indiv][roh_class] += length
f.close()

print('\t'.join(['pop','indiv','roh','len']))
for pop in dico :
    for indiv in dico[pop] :
        for roh_class in all_roh_classes :
            if roh_class not in dico[pop][indiv] :
                length = 0
            else :
                length = dico[pop][indiv][roh_class]
            new_line = [pop, indiv, roh_class, str(length)]
            print('\t'.join(new_line))
