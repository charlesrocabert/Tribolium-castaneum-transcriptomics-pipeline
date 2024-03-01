#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_ExtractHaplotypeBlocks.R
# --------------------------
# Extract haplotype-blocks from the genetic map.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Extract haplotype blocks from the genetic map ###
extract_haplotype_blocks <- function( map, margin_dist, min_size )
{
  HB_DATA = c()
  HB_name = 1
  for(CHR in unique(map$chr))
  {
    LG      = filter(map, chr==CHR)
    HB_size = table(LG$mean_pos)
    HB_size = sort(HB_size, decreasing=T)
    HB_size = HB_size[HB_size >= min_size]
    for(dist in as.numeric(names(HB_size)))
    {
      min_dist           = dist-margin_dist
      max_dist           = dist+margin_dist
      hb                 = filter(LG, mean_pos>=min_dist & mean_pos<=max_dist)
      hb$haplotype_block = rep(HB_name, dim(hb)[1])
      HB_DATA            = rbind(HB_DATA, hb)
      HB_name            = HB_name+1
    }
  }
  HB_DATA$ID = paste0(HB_DATA$chr,"-",HB_DATA$physical_pos)
  return(HB_DATA)
}


##################
#      MAIN      #
##################

#-----------------------------------------------#
# 1) Read command line arguments                #
#-----------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
SUFFIX          = args[4]
SIZE_THRESHOLD  = as.integer(args[5])
setwd(REPOSITORY_PATH)

#-----------------------------------------------#
# 2) Open the genetic map                       #
#-----------------------------------------------#
filename = paste0("./data/tribolium_ld/",POPULATION,"_",VERSION,"_",SUFFIX,".txt")
map      = read.table(filename, h=T, sep="\t")

#-----------------------------------------------#
# 3) Find haplotype blocks for each chromosome: #
#    ------------------------------------------ #
#    All SNPs with the same distance (in cM)    #
#    are considered to belong to the same       #
#    haplotype block.                           #
#    A minimum size is applied (e.g. minimum 10 #
#    SNPs per haplotype block).                 #
#-----------------------------------------------#
HB_DATA = extract_haplotype_blocks(map, MARGIN, SIZE_THRESHOLD)

#-----------------------------------------------#
# 4) Save the dataset                           #
#-----------------------------------------------#
filename_rds = paste0("./data/tribolium_haplotype_blocks/haplotype_blocks_",POPULATION,"_",VERSION,".rds")
filename_csv = paste0("./data/tribolium_haplotype_blocks/haplotype_blocks_",POPULATION,"_",VERSION,".csv")
saveRDS(HB_DATA, file=filename_rds)
write.table(HB_DATA, filename_csv, sep=";", col.names=T, row.names=F, quote=F)

