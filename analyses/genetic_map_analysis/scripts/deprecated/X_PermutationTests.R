#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_PermutationTests.R
# --------------------
# Run a permutation test on a given haplotype block.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Run SNP permutation test ###
run_snp_permutation_test <- function( dataset, variable, snp_list, nb_reps )
{
  x              = dataset[,variable]
  pos            = which(dataset$ID%in%snp_list)
  original_score = sum(x[pos])
  permutations   = t(sapply(1:nb_reps,function(i){
    score = sum(sample(x)[pos])
    return(score)
  }))
  qval = sum(permutations<original_score)/nb_reps
  return(1-qval)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#--------------------------------#
# 1) Load and prepare datasets   #
#--------------------------------#
DATA     = readRDS("./analyses/genetic_map_analysis/data/merged_dataset.rds")
HB_NAMES = unique(DATA$haplotype_block)
length(HB_NAMES)

sort(table(filter(DATA, evol_response_pool==1)$haplotype_block))

#--------------------------------#
# 3) Run the permutation test    #
#--------------------------------#
snps = filter(DATA, haplotype_block==559)$ID

run_snp_permutation_test(DATA, "isShifting", snps, 1000)

