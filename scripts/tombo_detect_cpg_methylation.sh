#!/bin/bash
#PBS -l nodes=1:ppn=33,mem=80gb
#PBS -l walltime=24:00:00

# Usage qsub -F "/fast5/dir stats_basename corr_group" \
# tombo_detect_cpg_methylation.sh

module load anaconda3
set -x
cd $PBS_O_WORKDIR

FAST5=$(realpath $1)
OUT_BASENAME=$2
CORRECTED_GROUP=$3

tombo detect_modifications alternative_model \
    --fast5-basedirs $FAST5 \
    --corrected-group $CORRECTED_GROUP \
    --alternate-bases CpG \
    --statistics-file-basename $OUT_BASENAME \
    --per-read-statistics-basename $OUT_BASENAME \
    --processes 32
