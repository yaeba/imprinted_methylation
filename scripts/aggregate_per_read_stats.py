#!/usr/bin/env python3
## Script to aggregate per-read statistics into coverage and proportion of reads
## supporting methylation call at the position
## Usage: python aggregate_per_read_stats.py [/path/per_read_df.pkl] out_prefix threshold(s)

import sys
import os
import numpy as np
import pandas as pd
import pickle

os.environ['R_HOME']='/stornext/System/data/apps/R/R-3.5.1/lib64/R'
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def open_pickled(file):
    with open(file, "rb") as fp:
        return pickle.load(fp)

def write_pickle(data, filename):
    with open(filename, "wb") as fp:
        pickle.dump(data, fp, protocol=pickle.HIGHEST_PROTOCOL)

def write_rdata(data, varname):
    r_data = pandas2ri.py2ri(data)
    robjects.r.assign(varname, r_data)
    robjects.r("save({}, file='{}')".format(varname, varname+".RData"))

i = 1
num_pkl = len(sys.argv) - 3
pickled = []

for _ in range(num_pkl):
    pickled.append(sys.argv[i])
    i += 1

out_prefix = sys.argv[i]
i += 1
thresh = eval(sys.argv[i])
try:
    len(thresh)
    single_threshold = False
except TypeError:
    single_threshold = True
    thresh = float(thresh)

print("Reading in pickled file(s)")
dfs = []
for p in pickled:
    dfs.append(open_pickled(p))
per_read_df = pd.concat(dfs, axis=0, ignore_index=True)
del dfs

print("Aggregating per read statistics")
per_read_df['coverage'] = 0

if single_threshold:
    # signal more likely to match against alternate base if modified LLR < thresh
    per_read_df['meth_prop'] = per_read_df['log_lik_ratio'].lt(thresh).astype(float)
else:
    per_read_df = per_read_df.loc[~per_read_df['log_lik_ratio'].between(thresh[0], thresh[1]),]
    per_read_df['meth_prop'] = per_read_df['log_lik_ratio'].lt(thresh[1]).astype(float)

aggregated_df = per_read_df.groupby(['chr', 'pos']).agg({'coverage': 'size',
                                                        'meth_prop': np.mean})
aggregated_df = aggregated_df.reset_index()

del per_read_df

print("Saving to pkl")
write_pickle(aggregated_df, out_prefix + ".pkl")
print("Saving to RData")
write_rdata(aggregated_df, out_prefix)
