#!/bin/bash
#PBS -l nodes=1:ppn=16,mem=50gb
#PBS -l walltime=72:00:00

# Usage: qsub -F "<mouse_per_read_stats.RData>" detect_dmr_dss.sh

set -x
cd $PBS_O_WORKDIR

module load R
export R_LIBS_USER="/stornext/Home/data/allstaff/t/tay.x/R/x86_64-pc-linux-gnu-library/3.5:$R_LIBS_USER"

RSCRIPT="/stornext/HPCScratch/home/tay.x/imprinted_methylation/scripts/detect_dmr_dss.R"
PER_READ_STATS=$(realpath $1)

Rscript $RSCRIPT $PER_READ_STATS
