#!/bin/bash

# Usage: bash nanopolish_phase_reads.sh reads.fastq mapped.bam genome.fa \
# vcf prefix

module load samtools

READS=$1
BAM=$2
GENOME=$3
VCF=$4
PREFIX=$5

nanopolish phase-reads -t 32 -r $READS -b $BAM -g $GENOME $VCF | samtools sort -T $PREFIX.tmp -o $PREFIX.phased.sorted.bam

samtools index $PREFIX.phased.sorted.bam
