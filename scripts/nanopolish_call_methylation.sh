#!/bin/bash
#PBS -l nodes=1:ppn=33,mem=10gb
#PBS -l walltime=90:00:00

# Usage: qsub -F "<fasta> <bam> <ref_genome> <output>" \
# nanopolish_call_methylation.sh

set -x
cd $PBS_O_WORKDIR

FASTA=$(realpath $1)
BAM=$(realpath $2)
REF=$(realpath $3)
OUTPUT=$(realpath $4)

nanopolish call-methylation \
    -r $FASTA \
    -b $BAM \
    -g $REF \
    -q "cpg" \
    -t 32 \
    > $OUTPUT
