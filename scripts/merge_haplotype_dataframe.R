#!/user/bin/Rscript

# Usage: Rscript merge_haplotype_dataframe.R <b6_base> <b6_signal> \
# <cast_base> <cast_signal> <output_prefix>

library(tidyverse)
library(reticulate)

args <- commandArgs(trailingOnly=TRUE)

loadRData <- function(infile) {
          # loads an RData file and returns it
          # assume that the variable saved in same name as the file
          load(infile)
          get(ls()[ls() == sub("^([^.]*).*", "\\1", basename(infile))])
}

b6_base <- loadRData(args[1])
b6_signal <- loadRData(args[2])
cast_base <- loadRData(args[3])
cast_signal <- loadRData(args[4])
output_basename <- args[5]

b6 <- bind_rows(list(b6_base, b6_signal))
cast <- bind_rows(list(cast_base, cast_signal))

b6 <- b6 %>%
  mutate(read_id=substring(read_id, first=1, last=36)) %>%
  group_by(read_id) %>%
  summarise(b6_score=mean(score, na.rm=TRUE),
  	    b6_num_align=n()) %>%
  ungroup()

cast <- cast %>%
  mutate(read_id=substring(read_id, first=1, last=36)) %>%
  group_by(read_id) %>%
  summarise(cast_score=mean(score, na.rm=TRUE),
            cast_num_align=n()) %>%
  ungroup()

# join haplotype results
assign(output_basename, inner_join(b6, cast, by="read_id"))

py_save_object(get(output_basename), file=paste0(output_basename, ".pkl"), 
               pickle="pickle")
save(list=output_basename, file=paste0(output_basename, ".RData"))