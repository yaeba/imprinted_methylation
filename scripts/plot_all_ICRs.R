#!/user/bin/Rscript

# Usage: Rscript plot_all_ICRs.R <mouse_per_read_stats.RData> \
# <cpg.tsv> <ICR_plots.tsv>

library(tidyverse)
library(gridExtra)

source("/stornext/HPCScratch/home/tay.x/imprinted_methylation/scripts/plot_methylation.R")

args <- commandArgs(trailingOnly=TRUE)

loadRData <- function(infile) {
          # loads an RData file and returns it
          # assume that the variable saved in same name as the file
          load(infile)
          get(ls()[ls() == sub("^([^.]*).*", "\\1", basename(infile))])
}

mouse_per_read_stats <- loadRData(args[1]) %>%
  dplyr::rename(start=pos) %>%
  mutate(end=start+1) %>%
  as.data.table()

cpg <- read_tsv(args[2], col_types="cii", col_names=FALSE)
names(cpg) <- c("chr", "start", "end")
cpg <- cpg %>%
  mutate(chr=sub("^chr", "", chr))

icr <- read_tsv(args[3], col_types="ciiccc")

# extract b6xcast and castxb6
mouse_per_read_stats <- mouse_per_read_stats %>%
  mutate(is_b6xcast=grepl('b6_maternal|cast_paternal', mouse_per_read_stats$strain),
         strain=sub("b6_|cast_", "", strain))
b6xcast <- mouse_per_read_stats %>%
  filter(is_b6xcast)
castxb6 <- mouse_per_read_stats %>%
  filter(!is_b6xcast)

layout <- rbind(c(1, 2), c(1, 2), c(3, 3), c(4, 4))

for (i in 1:nrow(icr)) {
  row <- icr[i, ]
  b6xcast_title <- paste0("Black6xCast_", row$name,
                                     ", imprints on ", row$allele,
                                     "aternal allele")
  castxb6_title <- paste0("CastxBlack6_", row$name,
                                     ", imprints on ", row$allele,
                                     "aternal allele")

  pr1 <- plotReadMethylation(b6xcast, row$chr, row$start, row$end,
     			    title=b6xcast_title)
  pr2 <- plotReadMethylation(castxb6, row$chr, row$start, row$end,
     			    title=castxb6_title)
  p1 <- plotMethylation(b6xcast, cpg, row$chr, row$start, row$end,
                        title=b6xcast_title)
  p2 <- plotMethylation(castxb6, cpg, row$chr, row$start, row$end,
                        title=castxb6_title)
  p <- grid.arrange(pr1, pr2, p1, p2, layout_matrix=layout)
  ggsave(paste0("ICR_", sub('/', '_', row$name), ".pdf"), plot=p, width=10, height=10)
  dev.off()
}