#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# plot_LG_SNP_repartition.R
# -------------------------
# Plot the repartition of SNPs among linkage groups.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")

### Plot the linkage groups ###
plot_linkage_groups <- function( data )
{
  data           = data[order(data$LG),]
  single_p       = sum(data$LG==0)/dim(data)[1]
  chrom_p        = sum(data$LG %in% 1:10)/dim(data)[1]
  text_label     = paste0(round(single_p,5), "% single SNPs,\n",round(chrom_p,5), "% SNPs in 10 first LGs")
  data_no_single = filter(data, LG != 0)
  p1 = ggplot(filter(data_no_single, LG < 21), aes(x=as.factor(LG))) +
    geom_bar(aes(y=(..count..)/sum(..count..))) +
    annotate("text", x=13, y=0.1, label=text_label, hjust=0, size=4) +
    ggtitle("SNPs repartition in first 20 LGs") +
    xlab("LG") +
    ylab("Frequency") +
    theme_classic()
  p2 = ggplot(filter(data_no_single, LG < 11), aes(x=as.factor(LG), fill=CHR)) +
    geom_bar(position="stack") +
    ggtitle("Chrom repartition in the 10 LGs") +
    xlab("LG") +
    ylab("Count") +
    theme_classic()
  theme_classic()
  p = plot_grid(p1, p2, ncol=2)
  return(p)
  #return(p2)
}


##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
setwd(REPOSITORY_PATH)

#--------------------------------#
# 2) Build the figures           #
#--------------------------------#

#LODlimit = 12
d = read.table(paste0("./data/tribolium_ld/stage_4_map.txt"), h=T, sep="\t")

p1 = plot_linkage_groups(d)
res = as.data.frame(d[!duplicated(d$CHR_LEN),])
res$COUNT = as.numeric(table(d$CHR))

p2 = ggplot(res, aes(CHR_LEN, COUNT, color=CHR)) +
  geom_point() +
  ggtitle("Chromosome length VS. Nb of SNPs per linkage group") +
  xlab("Chromosome length") +
  ylab("Nb. Snps per LG") +
  theme_classic() +
  theme(legend.position="none")

summary(lm(res$COUNT~res$CHR_LEN))
plot_grid(p1, p2, ncol=1)

###############################
###############################

# d1 = read.table(paste0("./data/tribolium_ld/LODlimit_exploration/LOD11_merged_map.txt"), h=T, sep="\t")
# d2 = read.table(paste0("./data/tribolium_ld/LODlimit_exploration/LOD12_merged_map.txt"), h=T, sep="\t")
# d3 = read.table(paste0("./data/tribolium_ld/LODlimit_exploration/LOD13_merged_map.txt"), h=T, sep="\t")
#
# p1 = plot_linkage_groups(d1)
# p2 = plot_linkage_groups(d2)
# p3 = plot_linkage_groups(d3)
#
# plot_grid(p1, p2, p3, ncol=3, labels="AUTO")







