#!/usr/bin/env python3

import sys

f = open('86_samples_ref_ancestral_LOF.recode.vcf')
positions = {}
for line in f :
    if line.startswith('#') :
        continue
    line = line.split()
    chrom, pos = int(line[0]), int(line[1])
    if chrom not in positions :
        positions[chrom] = []
    positions[chrom].append(pos)
f.close()


nb_block = 1000
genome_length = 0
for chrom in positions :
    start = positions[chrom][0]
    end = positions[chrom][-1]
    genome_length += end - start

print(f'Genome length: {genome_length}', file=sys.stderr)

def write_block(chrom, idx_block, block) :
    for pos in block :
        print(f'{idx_block}\t{chrom}\t{pos}')

block_size = genome_length / nb_block
block_size = 2399000            # hard coded to obtain 1000 blocks
chroms = list(positions.keys())
idx_block = 0
print('idx\tchrom\tpos')
for chrom in chroms :
    print(f'Working on chrom {chrom}', file=sys.stderr)
    idx_block += 1
    cur_block = [positions[chrom][0]]
    for idx in range(1, len(positions[chrom])) :
        pos = positions[chrom][idx]
        # append positions as long as block is not long enough
        if pos - cur_block[0] < block_size :
            cur_block.append(pos)
            continue
        # here, we need to write the current block and start a new one
        write_block(chrom, idx_block, cur_block)
        idx_block += 1
        cur_block = [pos]
    # here we still have a block to write (last of the chromosome)
    write_block(chrom, idx_block, cur_block)
    
    
print(block_size,file=sys.stderr)
print(idx_block,file=sys.stderr)
