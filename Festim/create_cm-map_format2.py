#!/usr/bin/env python

import gzip
import sys

for idx_chrom in range(1, 23) :
    print(f'Working on chrom {idx_chrom}', file=sys.stderr)
    f_in = open(f'genetic_maps/plink.chr{idx_chrom}.GRCh37.map')
    f_out = open(f'genetic_maps/map_chr{idx_chrom}.gmap', 'w')
    print('pposition rrate gposition', file=f_out)

    line = next(f_in)
    line = line.split()
    s_prev_recomb = line[2]
    prev_pos, prev_recomb = int(line[3]), float(line[2])
    print(f'{prev_pos} 0.0 {s_prev_recomb}', file=f_out)
    for line in f_in :
        line = line.split()
        s_cur_recomb = line[2]
        cur_pos, cur_recomb = int(line[3]), float(line[2])
        if cur_pos == prev_pos :
            continue
        if s_cur_recomb == s_prev_recomb :
            rrate = 0.0
        else :
            rrate = 1e6 * (cur_recomb - prev_recomb) / (cur_pos - prev_pos)
        if rrate < 0 :
            print(f'big error {rrate} {cur_pos}', file=sys.stderr)
            sys.exit(1)
        s_prev_recomb, prev_pos, prev_recomb = s_cur_recomb, cur_pos, cur_recomb
        print(f'{prev_pos} {rrate} {s_prev_recomb}', file=f_out)
    
    f_out.close()
    f_in.close()
    
