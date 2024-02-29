#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_GeneratePlinkFstClusters.R
# ----------------------------
# Generate the population clusters for PLINK 2.0.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation")

#-----------------------------------#
# 1) Loading sample list            #
#-----------------------------------#
samples = read.table("./data/tribolium_bam/samples_CT_HD_Tcas3.30.csv", h=T, sep=";")
LINES   = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2")

#-----------------------------------#
# 2) Generating population clusters #
#-----------------------------------#
clusters = c()
for(i in 1:dim(samples)[1])
{
  sample     = samples$sample[i]
  line       = samples$line[i]
  env        = samples$source_env[i]
  generation = samples$generation[i]
  if (line%in%LINES)
  {
    line     = c(sample, paste0(line,"-",env,"-",generation))
    clusters = rbind(clusters, line)
  }
}
clusters = as.data.frame(clusters)
names(clusters) = c("IID", "pop")
write.table(clusters, file="./data/tribolium_diversity/fst/CT_HD_Tcas3.30_clusters.txt", col.names=T, row.names=F, quote=F)

