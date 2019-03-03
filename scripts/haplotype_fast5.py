#!/usr/bin/env python3
# Script to split reads to corresponding haplotyped strain
# Usage python haplotype_fast5.py haplotype_df.pkl path/to/all_fast5.txt

import os
import sys
import h5py
import pickle
import tqdm
import re

def extract_read_id(fast5_file):
    """Function to extract read id from hdf5 file
    from a fast5 file."""
    assert os.path.isfile(fast5_file), "{} not found".format(fast5_file)
    with h5py.File(fast5_file, 'r') as fp:
        fastq = fp['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq']
        return fastq[()][1:50].decode('utf-8').split('_')[0]

pickled = sys.argv[1]
all_fast5 = sys.argv[2]

print("Reading pickled file")
with open(pickled, "rb") as file:
    haplotype_df = pickle.load(file)
    
print("Processing haplotype_df")
pattern = re.compile(".*score")
strains = list(filter(pattern.match, haplotype_df.columns))
strains = [x[:-6] for x in strains]

for strain in strains:
    os.system("mkdir -p {}".format(strain))

haplotype_df.set_index('read_id', inplace=True)


def final_genotype(diff):
    strain = None
    if diff > 0.2:
        strain = strains[0]
    elif diff < -0.2:
        strain = strains[1]
    return strain

haplotype_df['score_diff'] = haplotype_df[strains[0]+'_score'] - haplotype_df[strains[1]+'_score']
haplotype_df['final_genotype'] = haplotype_df['score_diff'].map(final_genotype)

genotype_dict = dict(zip(haplotype_df.index, haplotype_df['final_genotype']))

del haplotype_df

print("Processing all reads")


with tqdm.tqdm(total=os.path.getsize(all_fast5)) as pbar:
    with open(all_fast5, 'r') as fp:
        for line in fp:
            pbar.update(len(line))
            fast5 = line.rstrip()
            assert os.path.isfile(fast5), "{} not found".format(fast5)
        
            # to create folder to hold read (eg .../2/fast5)
            folder = fast5.split('/')[-2]
            read_id = extract_read_id(fast5)
            
            if read_id not in genotype_dict:
                continue
            
            strain = genotype_dict[read_id]
            
            path = "{}/{}/".format(strain, folder)
        
            if strain is not None:
                os.system("mkdir -p {}".format(path))
                os.system("cp {} {}/".format(fast5, path))
            
