#!/user/bin/Rscript

# Usage: Rscript bind_strains_dataframes.R <b6_mat> <cast_mat> <b6_pat> <cast_pat> \
# <output_basename>

library(tidyverse)

## Parse command
args <- commandArgs(trailingOnly=TRUE)

loadRData <- function(infile) {
          # loads an RData file and returns it
          # assume that the variable saved in same name as the file
          load(infile)
          get(ls()[ls() == gsub(basename(infile), pattern=".RData$", 
                                replacement="")])
}

b6_mat <- loadRData(args[1])
cast_mat <- loadRData(args[2])
b6_pat <- loadRData(args[3])
cast_pat <- loadRData(args[4])
output_basename <- args[5]

# append the strain information
b6_mat$strain <- "b6_maternal"
cast_mat$strain <- "cast_maternal"
b6_pat$strain <- "b6_paternal"
cast_pat$strain <- "cast_paternal"

# combine into a single dataframe
assign(output_basename, bind_rows(list(b6_mat, cast_mat, b6_pat, cast_pat)))

# compute the beta value (single-read single-site probabilty of methylation)
assign(output_basename, 
	get(output_basename) %>%
		mutate(beta=round(1/(1+exp(log_lik_ratio)), 6)))

# save the resultant dataframe
save(list=output_basename, file=paste0(output_basename, ".RData")) 