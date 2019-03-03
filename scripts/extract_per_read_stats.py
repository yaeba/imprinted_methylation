#!/usr/bin/env python3
## Script to load per-read stats file generated by tombo detect_modifications
## and make it into a dataframe with columns [chr, pos, log_lik_ratio, read_id]
## Usage: python extract_per_read_stats.py /path/per_read_stats \
## out_prefix <include_read_id>

import sys
import os
import numpy as np
import pandas as pd
import pickle
import itertools
from tombo import tombo_stats, tombo_helper

os.environ['R_HOME']='/stornext/System/data/apps/R/R-3.5.1/lib64/R'
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def write_pickle(data, filename):
    with open(filename, "wb") as fp:
        pickle.dump(data, fp, protocol=pickle.HIGHEST_PROTOCOL)

def write_rdata(data, varname):
    r_data = pandas2ri.py2ri(data)
    robjects.r.assign(varname, r_data)
    robjects.r("save({}, file='{}')".format(varname, varname+".RData"))

assert(len(sys.argv) >= 3)
infile = sys.argv[1]
out_prefix = sys.argv[2]
include_read_id = True if len(sys.argv) > 3 else False

print("Reading in per read stats")
# read in per-read statistics file
stats = tombo_stats.PerReadStats(infile)

print("Gathering statistics from all blocks")

# iterate and gather all statistics at CpG sites
rows = []

for chrm, strand, block_start, _, block_stats in stats:
    block_name = stats.blocks_index[(chrm, strand)][block_start]
    if len(block_stats) == 0:
        continue

    pos, log_lik_ratio, read_ids = np.array(list(zip(*block_stats)))
    pos += 1
    log_lik_ratio = np.round(log_lik_ratio, 3)
    if strand == '-':
        pos -=1
    
    if not include_read_id:
        rows.extend(list(zip(itertools.repeat(chrm), pos, log_lik_ratio)))
    else:
        # need to search for the read id
        if 'read_id_vals' in stats.per_read_blocks[block_name]:
            block_read_id_lookup = dict([
                (read_id_val, read_id) for read_id, read_id_val in
                zip(stats.per_read_blocks[block_name]['read_ids'].value,
                    stats.per_read_blocks[block_name]['read_id_vals'].value)])
        else:
            # read_ids previously stored (inefficiently as attributes)              
            block_read_id_lookup = dict([
                (read_id_val, read_id) for read_id, read_id_val in
                stats.per_read_blocks[block_name]['read_ids'].attrs.items()])

        read_id_vals = [eval(block_read_id_lookup[read_id]).decode('utf-8')
                        for read_id in read_ids]

        rows.extend(list(zip(itertools.repeat(chrm), pos, log_lik_ratio, read_id_vals)))

print("Storing to dataframe")

# create a dataframe
if not include_read_id:
    df = pd.DataFrame(rows, columns=['chr', 'pos', 'log_lik_ratio'])
else:
    df = pd.DataFrame(rows, columns=['chr', 'pos', 'log_lik_ratio', 'read_id'])

df = df.sort_values(['chr', 'pos'])
df['pos'] = df['pos'].astype(int)
df = df.reset_index(drop=True)

# save the dataframe in pkl and RData format
print("Saving to pkl")
write_pickle(df, out_prefix + ".pkl")
print("Saving to RData")
write_rdata(df, out_prefix)