#!/bin/bash
# Script to submit tombo runs to qsub

WRAPPER="/stornext/Bioinf/data/lab_speed/nanopore/keniry_mouse/txk/tombo/wrapper_no_filter.sh"
GENOME_B6FVB="/stornext/Bioinf/data/lab_speed/nanopore/keniry_mouse/genome_data/B6FVB_NJ.mgp.v5.snps.dbSNP142.fa"
GENOME_CAST="/stornext/Bioinf/data/lab_speed/nanopore/keniry_mouse/genome_data/CAST_EiJ.mgp.v5.snps.dbSNP142.fa"
B6FVB_READS="/stornext/HPCScratch/home/tay.x/Black6xCast_fast5_pass/B6FVB"
CAST_READS="/stornext/HPCScratch/home/tay.x/Black6xCast_fast5_pass/CAST"

SCRIPT="/stornext/HPCScratch/home/tay.x/scripts/run_tombo_mouse.sh"

cd /stornext/HPCScratch/home/tay.x/Black6xCast_fast5_pass

# run detect_modifications (alternative model and de novo) in parallel
qsub -F "$GENOME_B6FVB $B6FVB_READS B6FVB B6FVB_NJ 0" $SCRIPT

qsub -F "$GENOME_B6FVB $B6FVB_READS B6FVB B6FVB_NJ 1" $SCRIPT

qsub -F "$GENOME_CAST $CAST_READS CAST CAST_EiJ 0" $SCRIPT

qsub -F "$GENOME_CAST $CAST_READS CAST CAST_EiJ 1" $SCRIPT