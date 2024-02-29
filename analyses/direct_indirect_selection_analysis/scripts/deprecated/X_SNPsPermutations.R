#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_SNPsPermutations.R
# --------------------
# Parallelism permutation tests on SNPs.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("cowplot")

### Load eQTLs ###
load_eQTLs <- function( phenotype )
{
  #--------------------------#
  # 1) Load eQTLs data       #
  #--------------------------#
  EQTLs      = readRDS(paste0("./data/tribolium_eqtl/significant/HD_G1_Tcas3.30_imputed_",phenotype,"_significant.rds"))
  EQTLs$ID   = EQTLs$rs
  EQTLs$chr2 = as.numeric(as.factor(EQTLs$chr))
  #--------------------------#
  # 2) Load SNP annotation   #
  #--------------------------#
  ANNOTATION        = read.table("./data/tribolium_snp/snp_table_ALL_Tcas3.30_raw_SNP.csv", h=T, sep="\t")
  ANNOTATION$gene   = ANNOTATION$Feature_id
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "Annotation", "Putative_impact", "gene"))
  #--------------------------#
  # 3) Select and merge data #
  #--------------------------#
  eQTLs = merge(EQTLs, ANNOTATION_SELECT, by="ID")
  #--------------------------#
  # 4) Return data           #
  #--------------------------#
  rm(ANNOTATION)
  return(eQTLs)
}

### Return the stars corresponding to pvalues ###
get_pvalue_stars <- function( pvalues )
{
  N     = length(pvalues)
  stars = rep("", N)
  for(i in 1:N)
  {
    if (pvalues[i] < 0.0001)
    {
      stars[i] = "****"
    }
    if (pvalues[i] >= 0.0001 & pvalues[i] < 0.001)
    {
      stars[i] = "***"
    }
    if (pvalues[i] >= 0.001 & pvalues[i] < 0.01)
    {
      stars[i] = "**"
    }
    if (pvalues[i] >= 0.01 & pvalues[i] < 0.05)
    {
      stars[i] = "*"
    }
    if (pvalues[i] >= 0.05 & pvalues[i] < 0.1)
    {
      stars[i] = "."
    }
    if (pvalues[i] >= 0.1)
    {
      stars[i] = "ns"
    }
  }
  return(paste(stars, collapse = " | "))
}

### Run SNP permutations ###
run_snp_permutation_tests <- function( dataset, variable, snp_list, nb_reps )
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

#------------------------------#
# 1) Load the data             #
#------------------------------#
SNP_dataset = readRDS("./analysis/direct_indirect_selection_analysis/data/SNP_dataset.rds")
eQTLs       = load_eQTLs("EXPRESSION")

#------------------------------#
# 2) Extract list of SNPs      #
#------------------------------#
all_genes = unique(SNP_dataset$gene)
############
significant_sg_pool       = unique(filter(SNP_dataset, significant_sg_pool==1)$gene)
significant_sg_pool_cis   = unique(filter(SNP_dataset, significant_sg_pool==1)$ID)
significant_sg_pool_trans = filter(eQTLs, phenotype%in%significant_sg_pool)$ID
############
highest_sg_pool       = unique(filter(SNP_dataset, highest_sg_pool==1)$gene)
highest_sg_pool_cis   = unique(filter(SNP_dataset, highest_sg_pool==1)$ID)
highest_sg_pool_trans = filter(eQTLs, phenotype%in%highest_sg_pool)$ID
############
evol_response_pool       = unique(filter(SNP_dataset, evol_response_pool==1)$gene)
evol_response_pool_cis   = unique(filter(SNP_dataset, evol_response_pool==1)$ID)
evol_response_pool_trans = filter(eQTLs, phenotype%in%evol_response_pool)$ID
############
common_pool       = unique(filter(SNP_dataset, in_common_pool==1)$gene)
common_pool_cis   = unique(filter(SNP_dataset, in_common_pool==1)$ID)
common_pool_trans = filter(eQTLs, phenotype%in%common_pool)$ID
############
hub_genes       = unique(filter(SNP_dataset, hub_gene==1)$gene)
hub_genes_cis   = unique(filter(SNP_dataset, hub_gene==1)$ID)
hub_genes_trans = filter(eQTLs, phenotype%in%hub_genes)$ID
############
SNP_dataset$is_eQTL   = as.numeric(SNP_dataset$gene%in%eQTLs$gene)
SNP_dataset$has_eQTLs = as.numeric(SNP_dataset$gene%in%eQTLs$phenotype)

#------------------------------#
# 3) Run SNP permutation tests #
#------------------------------#
pval1 = run_snp_permutation_tests(SNP_dataset, "isFullyParallel", significant_sg_pool_cis, 10000)
pval2 = run_snp_permutation_tests(SNP_dataset, "isFullyParallel", highest_sg_pool_cis, 10000)
pval3 = run_snp_permutation_tests(SNP_dataset, "isFullyParallel", evol_response_pool_cis, 10000)
pval4 = run_snp_permutation_tests(SNP_dataset, "isFullyParallel", common_pool_cis, 10000)
pval5 = run_snp_permutation_tests(SNP_dataset, "isFullyParallel", hub_genes_cis, 10000)

#------------------------------#
# 4) Save the results          #
#------------------------------#
D = data.frame(
  c("significant_sg_pool", "highest_sg_pool", "evol_response_pool", "common_pool", "hub_genes"),
  c(pval1, pval2, pval3, pval4, pval5)
)
names(D) = c("snp_pool", "pvalue")
write.table(D, "./analysis/direct_indirect_selection_analysis/data/SNP_parallelism_permutations.txt", row.names=F, col.names=T, quote=F)

