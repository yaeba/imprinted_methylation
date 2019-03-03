#!/user/bin/Rscript

# Usage: Rscript generate_dss_paired.R b6.mat cast.mat b6.pat cast.pat

library(tidyverse)
library(DSS)
library(bsseq)

args <- commandArgs(trailingOnly=TRUE)

loadRData <- function(infile) {
	  # loads an RData file and returns it
	  # assume that the variable saved in same name as the file
	  load(infile)
	  get(ls()[ls() == sub("^([^.]*).*", "\\1", basename(infile))])
}

b6.mat <- loadRData(args[1])
cast.mat <- loadRData(args[2])
b6.pat <- loadRData(args[3])
cast.pat <- loadRData(args[4])

cols <- c("chr", "pos", "N", "X")
names(b6.mat) <- cols
names(cast.mat) <- cols
names(b6.pat) <- cols
names(cast.pat) <- cols

# create a BSseq object
BSobj = makeBSseqData(list(b6.mat, b6.pat, cast.mat, cast.pat),
      		c("b6_maternal", "b6_paternal", "cast_maternal", "cast_paternal"))


# compute dml
dml_parent <- DMLtest(BSobj, group1=c("b6_maternal", "cast_maternal"),
	      group2=c("b6_paternal", "cast_paternal"),
	      equal.disp=TRUE, smoothing=TRUE)

dml_strain <- DMLtest(BSobj, group1=c("b6_maternal", "cast_paternal"),
	      group2=c("cast_maternal", "b6_paternal"),
	      equal.disp=TRUE, smoothing=TRUE)

# compute dmr
dmr_parent <- callDMR(dml_parent, p.threshold=1e-4, dis.merge=1500)
dmr_strain <- callDMR(dml_strain, p.threshold=1e-4, dis.merge=1500)

# save the results
save(dml_parent, file="dml_parent.RData")
save(dml_strain, file="dml_strain.RData")
save(dmr_parent, file="dmr_parent.RData")
save(dmr_strain, file="dmr_strain.RData")
save(BSobj, file="BSobj.rda")
