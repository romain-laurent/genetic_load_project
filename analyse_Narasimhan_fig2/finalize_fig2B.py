#!/usr/bin/env python3

import sys
import random

def get_data(neutral, deleterious, pops) :
    print('Reading data...', file=sys.stderr)
    if 'CADD' in neutral :
        filename = 'counts_fig2B_CADD.txt'
    else :
        filename = 'counts_fig2B.txt'
    f = open(filename)
    next(f)                     # header
    dico = {}
    for line in f :
        chrom, pos, annot, pop, hom, alt = line.split()
        if annot != neutral and annot != deleterious :
            continue
        if (chrom, pos) not in dico :
            dico[(chrom, pos)] = {'annot':annot}
        dico[(chrom,pos)][pop] = (int(hom), int(alt))
        
    f.close()
    return dico

def count_hom_deleterious(data, deleterious, pops) :
    nb_hom = {}
    nb_alt = {}
    for pop in pops :
        nb_hom[pop] = 0
        nb_alt[pop] = {}
    for coords in data :
        annot = data[coords]['annot']
        if annot != deleterious :
            continue
        for pop in pops :
            tmp_hom, tmp_alt = data[coords][pop]
            if tmp_hom > 0 :
                nb_hom[pop] += 1
            if tmp_alt != 0 :
                if tmp_alt not in nb_alt[pop] :
                    nb_alt[pop][tmp_alt] = 0
                nb_alt[pop][tmp_alt] += 1
    return nb_hom, nb_alt

def prepare_neutral_data(data, pops, neutral) :
    dico = {}
    for pop in pops :
        dico[pop] = {}
    for coords in data :
        annot = data[coords]['annot']
        if annot != neutral :
            continue
        for pop in pops :
            tmp_hom, tmp_alt = data[coords][pop]
            if tmp_hom > 0 :
                tmp_hom = 1
            if tmp_alt != 0 :
                if tmp_alt not in dico[pop] :
                    dico[pop][tmp_alt] = []
                dico[pop][tmp_alt].append(tmp_hom)
    return dico

def perform_resampling(dico_hom, dico_freqs, pops) :
    result = {}
    for pop in pops :
        result[pop] = 0
        for freq in dico_freqs[pop] :
            count = dico_freqs[pop][freq]
            tmp = random.sample(dico_hom[pop][freq], count)
            result[pop] += sum(tmp)
    return result

def main() :
    pops = ['CTA','LPO','TUR','TJE']
    nb_resampling = 1000
    deleterious = 'CADD_sup_28'
    neutral = 'CADD_inf_10'
    data = get_data(neutral, deleterious, pops)
    nb_hom_deleterious, alt_freq_deleterious = count_hom_deleterious(data, deleterious, pops)
    dico_hom_neutral = prepare_neutral_data(data, pops, neutral)
    new_line = ['idx'] + pops
    print('\t'.join(new_line))
    new_line = ['0']
    for pop in pops :
        new_line.append(str(nb_hom_deleterious[pop]))
    print('\t'.join(new_line))
    for idx in range(nb_resampling) :
        new_line = [f'{idx+1}']
        resampling = perform_resampling(dico_hom_neutral, alt_freq_deleterious, pops)
        for pop in pops :
            new_line.append(str(resampling[pop]))
        print('\t'.join(new_line))
    
    
if __name__ == '__main__' :
    main()
