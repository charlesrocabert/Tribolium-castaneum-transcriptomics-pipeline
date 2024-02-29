#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_MergeAFs.R
# ------------
# Merge allelic frequencies for a given environment.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Collect all SNP allele frequencies ###
merge_allele_frequencies <- function( environment )
{
  GENERATIONS = c(1, 21)
  LINES       = c()
  POOLED      = c()
  if (environment == "CT")
  {
    LINES  = c("L1", "L3", "L5", "L6", "Mx1", "Mx3")
    POOLED = c("L4", "Mx1", "Mx3")
  } else if (environment == "HD")
  {
    LINES  = c("L1", "L2", "L3", "L5", "L6", "Mx1", "Mx2")
    POOLED = c("L2", "Mx1", "Mx2")
  }
  MERGED_DATA = c()
  for (generation in GENERATIONS)
  {
    for (line in LINES)
    {
      filename = paste0("./data/tribolium_snp/imputed/snp_table_", environment, "_G", generation, "_", line, "_Tcas3.30_imputed.csv")
      if (generation==1 & line%in%POOLED)
      {
        filename = paste0("./data/tribolium_snp/imputed/snp_table_", environment, "_G1_Tcas3.30_imputed.csv")
      }
      SNP            = read.table(filename, h=T, sep="\t", check.names=F)
      SNP            = SNP[,c("ID", "POS", "CHROM", "AF", "AC", "NCALLED")]
      SNP$NCALLED    = 2*SNP$NCALLED
      SNP$LINE       = rep(line, dim(SNP)[1])
      SNP$GENERATION = rep(generation, dim(SNP)[1])
      MERGED_DATA    = rbind(MERGED_DATA, SNP)
    }
  }
  return(MERGED_DATA)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

for (env in c("CT", "HD"))
{
  AF = merge_allele_frequencies(env)
  saveRDS(AF, file=paste0("./data/tribolium_afc/AF_",env,".rds"))
  #write.table(AF, file=paste0("./data/tribolium_afc/AF_",env,".csv"), sep=";", row.names=F, col.names=T, quote=F)
}

