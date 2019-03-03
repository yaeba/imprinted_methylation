#!/usr/bin/env python3
# Usage: python generate_haplotype_dataframe.py <bamfile> <vcf> <ref_genome> <out_prefix> \
# [offset, default=0]

import os
from collections import defaultdict
import warnings
import pysam
import numpy as np
import pandas as pd
import tqdm
import pickle
from functools import reduce
import sys

DEBUG = False

class HaplotypeDF(object):
    """ A class to create the haplotype dataframe """
    def __init__(self, bamfile, vcf, ref_genome):
        """Initialization
        Args:
            bamfile (str) : Path to the bamfile.
            vcf (str) : Path to the vcf file.
            ref_genome (str) : Path to the reference genome.
        """
        
        self.bamfile = pysam.AlignmentFile(bamfile, "rb")
        assert self.bamfile.check_index(), "Index not present for {}".format(bamfile)
        
        self.vcf = pysam.VariantFile(vcf)
        assert self.vcf.index is not None, "Index not present for {}".format(vcf)
        
        self.ref_genome = pysam.FastaFile(ref_genome)
        assert os.path.isfile(ref_genome+".fai"), "Index not present for {}".format(ref_genome)
        
        self.num_reads = self.count_num_reads(bamfile)

        
    def count_num_reads(self, bamfile):
        """Count the total number of reads mapped in bamfile"""
        num_reads = [int(stat.split('\t')[-2]) for stat in pysam.idxstats(bamfile).split('\n')[:-2]]
        return reduce(lambda x, y: x + y, num_reads)
      
    def haplotype_reads(self, invert=False, offset=0):
        """ Compute haplotyping scores for every read that mapped onto the reference genome
        
        Args:
            invert (bool) : Invert score to 1 - score (Used to prevent swapping columns of vcf).
            offset (float) : Offset for the basecall qualities in computing score.
        """

        with tqdm.tqdm(total=self.num_reads) as pbar:
            lst = []
            for aligned_read in self.bamfile.fetch():
                pbar.update(1)
                chrm = aligned_read.reference_name
                if not self.vcf.is_valid_reference_name(chrm):
                    # chrm/contig not in vcf, skip
                    continue

                if aligned_read.is_unmapped or aligned_read.is_secondary:
                        # use only primary and supplementary alignment
                        continue
                        
                # convert to 1-based for region
                start = aligned_read.reference_start + 1
                end = aligned_read.reference_end + 1

                # get a list of vcf overlapped in the aligned region
                snps = list(self.vcf.fetch(region="{}:{}:{}".format(chrm, start, end)))
                
                is_primary = 1 if not aligned_read.is_supplementary else 0

                # modify read query name to save space and take into account of multiple supp reads
                aligned_read.query_name = aligned_read.query_name[:36] + "_" + \
                                          ("primary" if is_primary else "supp_{}".format(aligned_read.reference_start))

                lst.append(np.array(self.assign_score(aligned_read, snps, invert, offset) + [is_primary]))

        self.haplotype_df = pd.DataFrame(lst, columns=['read_id', 'chrm', 'num_snps', 'ref_match', 
                                                   'alt_match', 'ref_score', 'alt_score', 'coverage',
                                                       'ratio', 'score', 'is_primary'])
        self.haplotype_df.dropna(axis=0, how='all', inplace=True)
        num_cols = ['num_snps', 'ref_match', 'alt_match', 'ref_score', 
                    'alt_score', 'coverage', 'ratio', 'score', 'is_primary']
        self.haplotype_df[num_cols] = self.haplotype_df[num_cols].apply(pd.to_numeric, errors="coerce")
        self.haplotype_df.reset_index(drop=True, inplace=True)
        
    def assign_score(self, read, snps, invert, offset):
        """
        Args:
            read (pysam.libcalignedsegment.AlignedSegment) : A single read.
            snps ([pysam.libcbcf.VariantRecord]) : List of SNPs that covered by read.
        """
  
        # ref and pos are 0-based
        ref_to_pos = {ref: pos for (pos, ref) in read.get_aligned_pairs() if ref is not None and pos is not None}
        
        num_snps = 0
        ref_match = 0
        alt_match = 0
        ref_score = 0
        alt_score = 0

        for snp in snps:
            # snp.pos is 1-based
            
            if (snp.pos-1) not in ref_to_pos.keys():
                # read does not cover snp
                continue

            if 'PASS' not in snp.filter.keys():
                # snp does not pass all the filters
                continue

            ref_base = self.ref_genome.fetch(region="{0}:{1}:{1}".format(read.reference_name, snp.pos))

            # indexing check
            if (not invert and ref_base != snp.ref) or (invert and ref_base not in snp.alts):
                print("Wrong indexing!")
                print("ref {} snp.alleles {} at pos {}".format(ref_base, snp.alleles, snp.pos))
                continue

            read_base_idx = ref_to_pos[snp.pos-1]
            read_base = read.query_sequence[read_base_idx]
            if read_base not in snp.alleles:
                # read something other than recognised SNP
                # this assumes ref_base is consistent with vcf file
                continue


            # with reference to paper by gigante
            quality = read.query_qualities[read_base_idx] + offset
            
            score = 1 - np.exp(-0.6927 - 0.1203 * quality)

            if read_base != snp.ref:
                alt_score += score
                alt_match += 1
            else:
                ref_score += score
                ref_match += 1
                
            num_snps += 1
        
        coverage = ref_score + alt_score
        
        if coverage == 0 or num_snps == 0:
            # read has no snp corresponding to vcf supplied
            return [read.query_name, read.reference_name, 0,
                             0, 0, 0, 0, 0, np.nan, np.nan]
    
        if invert:
            return [read.query_name, read.reference_name, num_snps,
                              alt_match, ref_match, round(alt_score, 2), 
                              round(ref_score, 2), round(coverage, 2), 
                              round(alt_score / coverage, 3),
                              round((alt_score + ref_match - ref_score) / num_snps, 3)]
        else:
            return [read.query_name, read.reference_name, num_snps,
                              ref_match, alt_match, round(ref_score, 2), 
                              round(alt_score, 2), round(coverage, 2), 
                              round(ref_score / coverage, 3),
                              round((ref_score + alt_match - alt_score) / num_snps, 3)]
     
    def save_pkl(self, filename="haplotype_df.pkl"):
        with open(filename, 'wb') as fp:
            pickle.dump(self.haplotype_df, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    def save_rdata(self, filename, varname):
        # write pandas dataframe to an .RData file
        import rpy2
        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        
        os.environ['R_HOME']='/stornext/System/data/apps/R/R-3.5.1/lib64/R'
        r_data = pandas2ri.py2ri(self.haplotype_df)
        robjects.r.assign(varname, r_data)
        robjects.r("save({}, file='{}')".format(varname, filename))
        

bamfile = sys.argv[1]
vcf = sys.argv[2]
ref_genome = sys.argv[3]
file_prefix = sys.argv[4]

offset = 0
if (len(sys.argv) > 5):
    offset = int(sys.argv[5])

df = HaplotypeDF(bamfile, vcf, ref_genome)
print("Computing scores for all reads")
print("Offset =", offset)
df.haplotype_reads(offset=offset)
df.save_pkl("{}.pkl".format(file_prefix))
df.save_rdata("{}.RData".format(file_prefix), file_prefix)
