#!/bin/bash
# Usage bash map.sh /path/to/genome.mmi reads output_prefix

module load minimap2
module load samtools

GENOME=$1
READS=$2
PREFIX=$3

# run minimap2 and sort
minimap2 -ax map-ont $GENOME $READS -t 32 --cs=long | samtools sort -T $PREFIX.tmp -o $PREFIX.sorted.bam

# index
samtools index $PREFIX.sorted.bam

