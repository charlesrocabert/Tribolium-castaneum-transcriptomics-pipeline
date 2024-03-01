#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# LODlimit_exploration_analysis.R
# -------------------------------
# Analyze the results of the LOD limit exploration.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")

### Get SNP proportions ###
get_SNP_proportions <- function( data )
{
  data     = data[order(data$LG),]
  single_p = sum(data$LG==0)/dim(data)[1]
  chrom_p  = sum(data$LG %in% 1:10)/dim(data)[1]
  sd_p     = sd(table(data$LG)[2:11]/sum(table(data$LG)[2:11]))
  res      = cbind(single_p, chrom_p, sd_p)
  colnames(res) = c("single", "chrom", "sd")
  return(res)
}

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
  p = plot_grid(p1, p2, ncol=1)
  return(p)
}


##################
#      MAIN      #
##################

#------------------------------------------#
# 1) Read command line arguments           #
#------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
setwd(REPOSITORY_PATH)

#------------------------------------------#
# 2) Plot the lodLimit exploration results #
#------------------------------------------#
LOD  = 5:20
prop = c()
for (lod in LOD)
{
  d = read.table(paste0("LODlimit_exploration/LOD",lod,"_merged_map.txt"), h=T, sep="\t")
  prop = rbind(prop, cbind(lod, get_SNP_proportions(d)))
}
par(mfrow=c(3,1))
plot(prop[,1], prop[,2], type="l", xlab="lodLimit", ylab="% SNPs", main="Proportion of single SNPs")
points(prop[,1], prop[,2], pch=20)
abline(v=10, col="red", lty=2)
plot(prop[,1], prop[,3], type="l", xlab="lodLimit", ylab="% SNPs", main="Proportion of SNPs in the 10 first LGs")
points(prop[,1], prop[,3], pch=20)
abline(v=10, col="red", lty=2)
plot(prop[,1], prop[,4], type="l", xlab="lodLimit", ylab="SD of relative proportions", main="Repartition heterogeneity in the 10 first LGs")
points(prop[,1], prop[,4], pch=20)
abline(v=10, col="red", lty=2)

#------------------------------------------#
# 3) Plot lodLimit=10 results              #
#------------------------------------------#
d1 = read.table("LOD10_map.txt", h=T, sep="\t")
d2 = read.table("LOD10_map_joined.txt", h=T, sep="\t")
d3 = read.table("LOD10_map_edited.txt", h=T, sep="\t")
p1 = plot_linkage_groups(d1) + ggtitle("Raw map")
p2 = plot_linkage_groups(d2) + ggtitle("Joined single SNPs map")
p3 = plot_linkage_groups(d3) + ggtitle("Cleaned map")
plot_grid(p1, p2, p3, ncol=3)

#------------------------------------------#
# 4) Comparison with chromosomes size      #
#------------------------------------------#
d    = read.table("LOD10_map_edited.txt", h=T, sep="\t")
info = d[!duplicated(d$LG),]
info = info[order(info$LG),]
RES  = table(d$LG)
RES  = cbind(RES, info$CHR_LEN)
plot(RES[,1], RES[,2], pch=20, xlab="Nb SNPs per LG", ylab="Chrom length")

