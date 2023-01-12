#!/usr/bin/env python3

pops = ['CTA','LPO','TUR','TJE']

authorized_gap = 0


def print_rohs(pop, indiv, all_rohs) :
    for roh in all_rohs :
        chrom, start, stop, cM = roh
        print(f'{pop} {indiv} {chrom} {start} {stop} {cM}')

cur_indiv = 'NA'
print('pop indiv chrom start stop cM')
for pop in pops :
    f = open(f'results_force_winsize/ROH_{pop}_winsize_50.roh.bed')
    for line in f :
        # when we start a new indiv
        if line.startswith('track') :
            if cur_indiv != 'NA' : # if we already have seen a real ROH
                all_rohs.append(prev_roh) # we record it
                # we print ROHs
                print_rohs(pop, cur_indiv, all_rohs)
            # in all cases, we start a new record
            line = line.split()
            cur_indiv = line[2]
            all_rohs = []
            prev_roh = [0, 0, 0, 0]
            continue
        line = line.split()
        chrom, start, end, cM = line[0][3:], int(line[1]), int(line[2]), float(line[4])
        # when we start a new chromosome
        if chrom != prev_roh[0] :
            if prev_roh[0] != 0 : # if we already have seen a real ROH
                all_rohs.append(prev_roh) # we record it as such (it cannot extend)
            # in all cases, we start a new record
            prev_roh = [chrom, start, end, cM]
            continue

        # in all other cases
        # if we need to extend, we do
        if start - prev_roh[2] <= authorized_gap :
            prev_roh[2] = end
            prev_roh[3] += cM
            continue
        # else, the previous ROH cannot extend -> we record it
        all_rohs.append(prev_roh)
        prev_roh = [chrom, start, end, cM]

    f.close()
    # we print the last indiv of the pop
    all_rohs.append(prev_roh) # we record it
    print_rohs(pop, cur_indiv, all_rohs)
    cur_indiv = 'NA'
    all_rohs = []
    prev_roh = [0, 0, 0, 0]
            
