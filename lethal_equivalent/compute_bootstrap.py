#!/usr/bin/env python3

import sys
import random

def get_args() :
    if len(sys.argv) != 3 :
        print(f'usage: {sys.argv[0]} neutral non_neutral', file=sys.stderr)
        sys.exit(1)
    neutral, non_neutral = sys.argv[1], sys.argv[2]
    if 'CADD' not in neutral and 'VEP' not in neutral :
        print('Unknown mutation type. Exiting...', file=sys.stderr)
        sys.exit(1)
    if 'CADD' in neutral :
        if not 'CADD' in non_neutral :
            print('Incompatible neutral and non neutral. Exiting...', file=sys.stderr)
            sys.exit(1)
        block_filename = '../bootstrap_blocks_definition.txt'
        geno_filename = '../../annotation/refined_86_samples_ref_is_ancestral.vcf'
        annotation_filename = '../CADD_scores_computations/scores_with_classes.txt'
    else :
        if not 'VEP' in non_neutral :
            print('Incompatible neutral and non neutral. Exiting...', file=sys.stderr)
            sys.exit(1)
        block_filename = '../bootstrap_VEP/bootstrap_blocks_definition.txt'
        geno_filename = '../bootstrap_VEP/86_samples_ref_ancestral_LOF.recode.vcf'
        annotation_filename = '../VEP/VEP_annotations.txt'

    print('I will be working with:', file=sys.stderr)
    print(f'\tneutral: {neutral}', file=sys.stderr)
    print(f'\tnon neutral: {non_neutral}', file=sys.stderr)
    print(f'\tblock definitions: {block_filename}', file=sys.stderr)
    print(f'\tgenotypes: {geno_filename}', file=sys.stderr)
    print(f'\tannotations: {annotation_filename}', file=sys.stderr)
    return neutral, non_neutral, block_filename, geno_filename, annotation_filename

def get_samples(pops) :
    print('-> reading samples infos...', file=sys.stderr)
    dico = {}
    for pop in pops :
        f = open(f'../samples_{pop}.txt')
        lines = f.readlines()
        f.close()
        dico[pop] = [line.strip() for line in lines]
    return dico

def get_blocks(filename) :
    print('-> reading blocks infos...', file=sys.stderr)
    dico = {}
    dico2 = {}
    block_names = set()
    f = open(filename)
    next(f)                     # header
    for line in f :
        line = line.split()
        block_name, chrom, pos = line
        if (chrom, pos) in dico :
            print('Error: duplicated SNP', file=sys.stderr)
        dico[(chrom, pos)] = [block_name]
        if block_name not in dico2 :
            dico2[block_name] = set()
        dico2[block_name].add((chrom, pos))
        block_names.add(block_name)
    f.close()
    block_names = list(block_names) # cast to list because the random functions used later to not suport a set
    return (dico, dico2, block_names)

def get_annotations(filename, dico) :
    print('-> reading annotations infos...', file=sys.stderr)
    f = open(filename)
    next(f)
    for line in f :
        line = line.split()
        chrom, pos, annot = line[0], line[1], line[-1]
        if (chrom, pos) not in dico :
            print('Error: annotated SNP not in any block', file=sys.stderr)
        dico[(chrom, pos)].append(annot)
    f.close()

def count_one_line(dico_idxs, line) :
    counts = {}
    for pop, idxs in dico_idxs.items() :
        counts[pop] = [0,0]     # number of hom derived, number of derived allele
        for idx in idxs :
            geno = line[idx]
            if geno == '0/0' :  # homozygous ancestral -> nothing to count
                continue
            elif geno == '1/1' :  # homozygous derived
                counts[pop][0] += 1
                counts[pop][1] += 2
            else :              # heterozygous
                counts[pop][1] += 1
    return counts
    
def count_genotypes(dico_snps, dico_samples, filename, pops) :
    print('-> reading genotypes and counting stuff...', file=sys.stderr)
    f = open(filename)
    counts = {}
    for line in f :
        if line.startswith('##') : # comments
            continue
        if line.startswith('#') : # header -> get idxs of indivs
            line = line.split()
            dico_idxs = {}
            for pop, samples in dico_samples.items() :
                dico_idxs[pop] = []
                for sample in samples :
                    dico_idxs[pop].append(line.index(sample))
            continue
        line = line.split()
        # here we go
        chrom, pos = line[0], line[1]
        tmp_counts = count_one_line(dico_idxs, line)
        block_id = dico_snps[(chrom,pos)][0]
        if block_id not in counts :
            counts[block_id] = []
        to_append = [(chrom, pos)]
        for pop in pops :
            to_append.append(tmp_counts[pop])
        counts[block_id].append(to_append)
    f.close()

    return counts


def compute_final_stat(dico_samples, pops, counts, dico_snps, neutral, non_neutral, wanted_blocks) :
    to_return = []
    for idx_pop, pop in enumerate(pops) :
        sample_size = len(dico_samples[pop])
        # the possible number of allele counts is 2*sample_size+1
        hom_per_freq = {neutral:[0]*(2*sample_size + 1), non_neutral:[0]*(2*sample_size + 1)}
        site_per_freq = {neutral:[0]*(2*sample_size + 1), non_neutral:[0]*(2*sample_size + 1)}
        for block_id in wanted_blocks : # for each block
            for snp in counts[block_id] : # for each SNP in block
                coords = snp[0]
                snp_annot = dico_snps[coords][1] # find annotation of this SNP
                if snp_annot != neutral and snp_annot != non_neutral :
                    continue    # if not interesting annotation, skip
                # otherwise, update relevant counters
                nb_hom, nb_alt = snp[idx_pop+1]
                hom_per_freq[snp_annot][nb_alt] += nb_hom
                site_per_freq[snp_annot][nb_alt] += 1

        # now we have the whole genome, time to finalize
        # for each allelic frequency (except 0 and 2*sample_size)
        final_depletion = 0
        for nb_alt in range(1, 2*sample_size) :
            nb_hom_neutral = hom_per_freq[neutral][nb_alt]
            nb_site_neutral = site_per_freq[neutral][nb_alt]
            nb_hom_non_neutral = hom_per_freq[non_neutral][nb_alt]
            nb_site_non_neutral = site_per_freq[non_neutral][nb_alt]
            if nb_site_non_neutral == 0 or nb_site_neutral == 0 :
                # other solution -> add one to denominator
                #nb_site_non_neutral = max(nb_site_non_neutral, 1)
                #nb_site_neutral = max(nb_site_neutral, 1)
                continue        # prevent division by 0 -> skip those cases
        
            # compute depletion
            mean_neutral = nb_hom_neutral / nb_site_neutral
            mean_non_neutral = nb_hom_non_neutral / nb_site_non_neutral
            depletion = mean_neutral - mean_non_neutral
            # weighted sum like Narasimhan
            final_depletion += depletion * nb_alt / (2*sample_size)
        to_return.append(final_depletion)
    return to_return
    
def main() :
    pops = ['CTA','LPO','TUR','TJE']
    nb_bootstrap = 1000
    
    # get command line args
    neutral, non_neutral, block_filename, geno_filename, annotation_filename = get_args()
    # get sample ids for each pop
    dico_samples = get_samples(pops)
    # get block id for each SNP
    dico_snps, dico_blocks, block_names = get_blocks(block_filename)
    # get annotations
    get_annotations(annotation_filename, dico_snps)
    # count stuff that will be used for final stats
    counts = count_genotypes(dico_snps, dico_samples, geno_filename, pops)

    # now we have all the sumstats we need to finalize computations and perform bootstrap
    print('-> computing stat on real data...', file=sys.stderr)
    depletions = compute_final_stat(dico_samples, pops, counts, dico_snps, neutral, non_neutral, block_names)
    
    new_line = ['idx'] + pops
    print('\t'.join(new_line))
    new_line = ['0'] + [str(i) for i in depletions]
    print('\t'.join(new_line))
    print('-> computing bootstrap...', file=sys.stderr)
    for idx_boot in range(nb_bootstrap) :
        wanted = random.choices(block_names, k=len(block_names))
        depletions = compute_final_stat(dico_samples, pops, counts, dico_snps, neutral, non_neutral, wanted)
        new_line = [f'{idx_boot+1}'] + [str(i) for i in depletions]
        print('\t'.join(new_line))
    

if __name__ == '__main__' :
    main()
