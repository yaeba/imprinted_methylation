#!/bin/bash
#PBS -l nodes=1:ppn=17,mem=80gb
#PBS -l walltime=72:00:00

# Script to run tombo detect_modifications on haplotyped mouse reads
# Usage qsub -F "ref_genome reads prefix corrected_group de_novo" run_tombo_mouse.sh

module load anaconda3
module load R

GENOME=$1
READS=$2
PREFIX=$3
CORRECTED_GROUP=$4
DE_NOVO=$5

WRAPPER="/stornext/Bioinf/data/lab_speed/nanopore/keniry_mouse/txk/tombo/wrapper_no_filter.sh"

cd /stornext/HPCScratch/home/tay.x/Black6xCast_fast5_pass

# run detect_modifications
if [ $DE_NOVO -eq 0 ]; then
	bash $WRAPPER -g $GENOME -f $READS -a -ab CpG -o -p $PREFIX.alt.CpG -c $CORRECTED_GROUP -pr 16
else
	bash $WRAPPER -g $GENOME -f $READS -d -o -p $PREFIX.de.novo -c $CORRECTED_GROUP -pr 16 
fi

echo "All done"
