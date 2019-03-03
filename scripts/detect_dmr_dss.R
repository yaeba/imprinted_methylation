#!/user/bin/Rscript

# Usage: Rscript detect_dmr_dss.R <mouse_per_read.RData>

library(tidyverse)
library(DSS)
library(bsseq)

args <- commandArgs(trailingOnly=TRUE)

# function to load saved RData to a separate variable name
loadRData <- function(infile) {
          # loads an RData file and returns it
          # assume that the variable saved in same name as the file
          load(infile)
          get(ls()[ls() == sub("^([^.]*).*", "\\1", basename(infile))])
}


print("Reading RData")

per_read_stats <- loadRData(args[1])

print("Aggregating per-read stats")
aggregated_stats <- per_read_stats %>%
	group_by(strain, chr, pos) %>%
	summarise(N=n(),
		  X=sum(beta > 0.5)) %>%
	ungroup()


print("Creating BSseq object")
# create a BSseq object
BSobj = makeBSseqData(list(
      aggregated_stats %>% filter(strain=='b6_maternal') %>% select(-c(strain)), 
      aggregated_stats %>% filter(strain=='b6_paternal') %>% select(-c(strain)),
      aggregated_stats %>% filter(strain=='cast_maternal') %>% select(-c(strain)),
      aggregated_stats %>% filter(strain=='cast_paternal') %>% select(-c(strain))),
      		c("b6_maternal", "b6_paternal", "cast_maternal", "cast_paternal"))


print("Computing dml")
# compute dml
dml_parent <- DMLtest(BSobj, group1=c("b6_maternal", "cast_maternal"),
	      group2=c("b6_paternal", "cast_paternal"),
	      equal.disp=TRUE, smoothing=TRUE)

dml_strain <- DMLtest(BSobj, group1=c("b6_maternal", "cast_paternal"),
	      group2=c("cast_maternal", "b6_paternal"),
	      equal.disp=TRUE, smoothing=TRUE)

print("Computing dmr")
# compute dmr
dmr_parent <- callDMR(dml_parent, p.threshold=1e-4, dis.merge=1500)
dmr_strain <- callDMR(dml_strain, p.threshold=1e-4, dis.merge=1500)

print("Saving the results to RData")
# save the results
save(dml_parent, file="dml_parent.RData")
save(dml_strain, file="dml_strain.RData")
save(dmr_parent, file="dmr_parent.RData")
save(dmr_strain, file="dmr_strain.RData")
save(BSobj, file="BSobj.rda")
