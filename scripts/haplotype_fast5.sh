#!/bin/bash
# Usage: bash haplotype_fast5.sh /path/to/haplotype_df.pkl n_cores [/paths/to/fast5s]


module load parallel

SCRIPT="/stornext/HPCScratch/home/tay.x/imprinted_methylation/scripts/haplotype_fast5.py"
PICKLED=$(realpath $1); shift
N_CORES=$1; shift

find $@ -name "*.fast5" > all_fast5.txt

n_lines=$(echo "$(cat all_fast5.txt | wc -l) / $N_CORES + 1" | bc)

split -l $n_lines all_fast5.txt splitted.

find . -maxdepth 1 -name "splitted.*" | \
    parallel -I% --max-args 1 python $SCRIPT $PICKLED %

rm -f splitted.*
