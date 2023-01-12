#!/usr/bin/env python3

pops = ['CTA','LPO','TUR','TJE']
wanted_winsize = 50
wanted_class = 'C'

f_out = open(f'ROH_class_{wanted_class}_winsize_{wanted_winsize}.txt', 'w')
print('pop indiv chrom start stop', file=f_out)
for pop in pops :
    filename = f'results_force_winsize/ROH_{pop}_winsize_{wanted_winsize}.roh.bed'
    f = open(filename)
    for line in f :
        if line.startswith('track') :
            line = line.split()
            indiv = line[2]
            continue
        line = line.split()
        chrom = line[0][3:]
        start, stop = line[1], line[2]
        if line[3] != wanted_class :
            continue
        new_line = [pop, indiv, chrom, start, stop]
        print(' '.join(new_line), file=f_out)
    f.close()
f_out.close()

