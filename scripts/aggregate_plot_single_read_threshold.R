#!/user/bin/Rscript

# Usage: Rscript aggregate_plot_single_read_threshold.R \
# <all_per_read_stats.RData> <plot_basename> <threshold_1> [threshold_2]

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
load(args[1])
all_per_read <- get(ls()[ls() != "args"])
plot_basename <- args[2]
thresh_1 <-as.numeric(args[3])
thresh_2 <- NULL
if (length(args) > 3) {
   thresh_2 <- as.numeric(args[4])
}

print("Aggregating")
aggregated <- all_per_read %>%
	   group_by(chr, pos, status) %>%
	   summarise(cov=ifelse(is.null(thresh_2),
			n(),
			sum(p_val<thresh_1 | p_val>thresh_2)),
		     meth_prop=sum(p_val<thresh_1)/cov) %>%	
	   ungroup() %>%
	   na.omit()

text <- ifelse(is.null(thresh_2), as.character(thresh_1),
     				     paste0("(", thresh_1, ", ", thresh_2, ")"))

print("Plotting")
p <- ggplot(aggregated) +
  geom_density(aes(x=meth_prop, fill=status), alpha=0.5, color="white", size=0.01) +
  theme_bw() +
  ggtitle(paste0("Fraction supporting methylation per site, threshold=", text))

ggsave(paste0(plot_basename, "_", text, ".pdf"),
       plot=p, width=7, height=4)
