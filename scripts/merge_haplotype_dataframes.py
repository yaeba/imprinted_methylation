#!/usr/bin/env python3
## Script to merge the haplotype dataframes generated
## Usage: python merge_haplotype_dataframes.py /path/b6_base_df \
## /path/b6_signal_df /path/cast_base_df /path/cast_signal_df output_prefix

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

if len(sys.argv) != 6:
    print("Invalid number of args")
    exit(1)

b6_base_path = sys.argv[1]
b6_signal_path = sys.argv[2]
cast_base_path = sys.argv[3]
cast_signal_path = sys.argv[4]
out_prefix = sys.argv[5]

## read the pickled files
b6_base = open_pickled(b6_base_path).set_index("read_id")
b6_signal = open_pickled(b6_signal_path).set_index("read_id")
cast_base = open_pickled(cast_base_path).set_index("read_id")
cast_signal = open_pickled(cast_signal_path).set_index("read_id")

## add prefixes to the columns
b6_base.columns = "b6_base_" + b6_base.columns
b6_signal.columns = "b6_signal_" + b6_signal.columns
cast_base.columns = "cast_base_" + cast_base.columns
cast_signal.columns = "cast_signal_" + cast_signal.columns

## merge the dataframes
haplotype_df = b6_base.merge(b6_signal, how="inner", left_index=True, right_index=True) \
                      .merge(cast_base, how="inner", left_index=True, right_index=True) \
                      .merge(cast_signal, how="inner", left_index=True, right_index=True)
haplotype_df = haplotype_df.reset_index()

## delete unused dataframes to save space
del b6_base, b6_signal, cast_base, cast_signal

## compute base and signal genotype
def basic_genotype(x):
    haplotype_value = np.sign(x - 0.5)
    if haplotype_value == 1:
        return "ref"
    elif haplotype_value == -1:
        return "alt"
    else:
        return "NA"

haplotype_df['b6_base_genotype'] = haplotype_df['b6_base_score'].map(basic_genotype)
haplotype_df['b6_signal_genotype'] = haplotype_df['b6_signal_score'].map(basic_genotype)
haplotype_df['cast_base_genotype'] = haplotype_df['cast_base_score'].map(basic_genotype)
haplotype_df['cast_signal_genotype'] = haplotype_df['cast_signal_score'].map(basic_genotype)

## add "info" and "genotype" columns
b6_info = list()
cast_info = list()
b6_genotype = list()
cast_genotype = list()

for _, row in haplotype_df.iterrows():
    
    if row['b6_base_genotype'] == "NA" and row['b6_signal_genotype'] == "NA":
        b6_info.append("low_coverage")
        b6_genotype.append("NA")
    elif row['b6_base_num_snps'] < 5 and row['b6_signal_num_snps'] < 5:
        b6_info.append("low_coverage")
        b6_genotype.append("fail")
    elif row['b6_base_genotype'] == row['b6_signal_genotype']:
        b6_info.append("pass")
        b6_genotype.append(row['b6_base_genotype'])
    elif row['b6_base_num_snps'] > 3*row['b6_signal_num_snps']:
        b6_info.append("signal_low_coverage")
        b6_genotype.append(row['b6_base_genotype'])
    elif row['b6_signal_num_snps'] > 3*row['b6_base_num_snps']:
        b6_info.append("base_low_coverage")
        b6_genotype.append(row['b6_signal_genotype'])
    elif np.abs(row['b6_base_score'] - 0.5) > 3*np.abs(row['b6_signal_score'] - 0.5):
        b6_info.append("signal_uncertain")
        b6_genotype.append(row['b6_base_genotype'])
    elif np.abs(row['b6_signal_score'] - 0.5) > 3*np.abs(row['b6_base_score'] - 0.5):
        b6_info.append("base_uncertain")
        b6_genotype.append(row['b6_signal_genotype'])
    else:
        b6_info.append("fail")
        b6_genotype.append("fail")
        

    if row['cast_base_genotype'] == "NA" and row['cast_signal_genotype'] == "NA":
        cast_info.append("low_coverage")
        cast_genotype.append("NA")
    elif row['cast_base_num_snps'] < 5 and row['cast_signal_num_snps'] < 5:
        cast_info.append("low_coverage")
        cast_genotype.append("fail")
    elif row['cast_base_genotype'] == row['cast_signal_genotype']:
        cast_info.append("pass")
        cast_genotype.append(row['cast_base_genotype'])
    elif row['cast_base_num_snps'] > 3*row['cast_signal_num_snps']:
        cast_info.append("signal_low_coverage")
        cast_genotype.append(row['cast_base_genotype'])
    elif row['cast_signal_num_snps'] > 3*row['cast_base_num_snps']:
        cast_info.append("base_low_coverage")
        cast_genotype.append(row['cast_signal_genotype'])
    elif np.abs(row['cast_base_score'] - 0.5) > 3*np.abs(row['cast_signal_score'] - 0.5):
        cast_info.append("signal_uncertain")
        cast_genotype.append(row['cast_base_genotype'])
    elif np.abs(row['cast_signal_score'] - 0.5) > 3*np.abs(row['cast_base_score'] - 0.5):
        cast_info.append("base_uncertain")
        cast_genotype.append(row['cast_signal_genotype'])
    else:
        cast_info.append("fail")
        cast_genotype.append("fail")
        
haplotype_df['b6_genotype'] = b6_genotype
haplotype_df['b6_info'] = b6_info
haplotype_df['cast_genotype'] = cast_genotype
haplotype_df['cast_info'] = cast_info

## save the merged dataframes in pkl and RData format
write_pickle(haplotype_df, out_prefix + ".pkl")
write_rdata(haplotype_df, out_prefix)
