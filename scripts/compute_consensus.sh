#!/bin/bash
#PBS -l nodes=1:ppn=25,mem=10gb,walltime=72:00:00

# Usage qsub -F "draft.fa reads.fasta reads.bam prefix" compute_consensus.sh



module unload nanopolish
module load anaconda3
module load parallel

MAKERANGE="/stornext/System/data/apps/nanopolish/nanopolish-0.8.5/scripts/nanopolish_makerange.py"

set -x
cd $PBS_O_WORKDIR

DRAFT=$(realpath $1)
READS=$(realpath $2)
BAM=$(realpath $3)
PREFIX=$4

# split genome in 50kb blocks and run in parallel
python $MAKERANGE $DRAFT | parallel --results $PREFIX.results -P 3 \
    nanopolish variants --consensus \
    -o $PREFIX.polished.{1}.vcf \
    -w {1} \
    -r $READS \
    -b $BAM \
    -g $DRAFT \
    -t 8 \
    --min-candidate-frequency 0.1
