#!/bin/bash
#PBS -l nodes=1:ppn=33,mem=80gb
#PBS -l walltime=36:00:00

# Usage qsub -F "/fast5/dir /path/ref_genome.fa corrected_group" \
# tombo_resquiggle.sh


module load anaconda3
set -x
cd $PBS_O_WORKDIR

FAST5=$(realpath $1)
REF=$(realpath $2)
CORRECTED_GROUP=$3

tombo resquiggle $FAST5 $REF --corrected-group $CORRECTED_GROUP \
    --dna --q-score 5.0 \
    --signal-matching-score 2.0 \
    --num-most-common-errors 10 \
    --failed-reads-filename $CORRECTED_GROUP.failed_reads.txt \
    --processes 32 \
    --overwrite
