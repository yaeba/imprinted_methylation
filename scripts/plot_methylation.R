#!/user/bin/Rscript

library(tidyverse)
library(rtracklayer)
library(data.table)


plotMethylation <- function(df, cpg, chrm, start, end, overhang=20000, 
                            title="Title") {
  region_start <- start - overhang
  region_end <- end + overhang
  df <- df %>%
    filter(chr==chrm, start>region_start, end<region_end)
  
  p <- ggplot(df) + 
    geom_smooth(method="loess", 
                aes(x=start, y=beta, group=read_id, color=strain),
                span=2, se=FALSE) +
    ggplot2::ylim(0, 1) +
    geom_point(data=cpg %>% filter(chr==chrm,
                                   start > region_start,
                                   end < region_end),
               aes(x=start), y=0, pch='|') +
    geom_ribbon(data=data_frame(x=c(start, end)),
                aes(x), fill='red', color='transparent',
                alpha=0.2, ymin=0, ymax=1) +
    
    scale_color_manual(values=c('maternal'='#F8766D',
                                'paternal'='#00BFC4')) +
    theme_bw() +
    ggtitle(title) +
    xlab(paste0("chr ", chrm, " pos"))
  p
}

plotReadMethylation <- function(df, chrm, start, end, overhang=2000, 
                                title="Title") {
  region_start <- start - overhang
  region_end <- end + overhang
  df <- df %>%
    filter(chr==chrm, start>region_start, end<region_end)
  p <- ggplot(df, aes(x=start, group=read_id, y=1, color=beta)) + 
    geom_line(position=ggstance::position_dodgev(height=0.2), size=5) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5, space="Lab") +
    geom_vline(xintercept=c(start, end), linetype="dashed") + 
    facet_grid(rows=vars(strain)) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(x=paste0("chr ", chrm, " pos"), y="read", title=title)
  p
}

