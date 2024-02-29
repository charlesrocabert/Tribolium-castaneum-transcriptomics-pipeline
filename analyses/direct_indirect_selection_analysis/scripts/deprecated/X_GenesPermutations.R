#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 4_GenePermutations.R
# --------------------
# Permutation tests on genes.
# (LOCAL SCRIPT)
#***************************************************************************

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

### Run gene sum permutations ###
run_gene_permutation_tests <- function( dataset, variable, gene_list, nb_reps )
{
  x              = dataset[,variable]
  pos            = which(dataset$gene%in%gene_list)
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

#----------------------------#
# 1) Load the data           #
#----------------------------#
gene_dataset = readRDS("./analysis/direct_indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset = readRDS("./analysis/direct_indirect_selection_analysis/data/eQTL_dataset.rds")
eQTLs        = load_eQTLs("EXPRESSION")

#----------------------------#
# 2) Extract list of genes   #
#----------------------------#
all_genes             = unique(gene_dataset$gene)
significant_sg_pool   = unique(filter(gene_dataset, significant_sg_pool==1)$gene)
highest_sg_pool       = unique(filter(gene_dataset, highest_sg_pool==1)$gene)
evol_response_pool    = unique(filter(gene_dataset, evol_response_pool==1)$gene)
common_pool           = unique(filter(gene_dataset, in_common_pool==1)$gene)
hub_genes             = unique(filter(gene_dataset, hub_gene==1)$gene)
eQTL_carriers         = unique(filter(gene_dataset, eQTL_carrier==1)$gene)
eQTL_phenotypes       = unique(filter(gene_dataset, eQTL_phenotype==1)$gene)
no_pleiotropy_genes   = unique(filter(eQTL_dataset, Pleiotropy_category=="No pleiotropy")$gene)
low_pleiotropy_genes  = unique(filter(eQTL_dataset, Pleiotropy_category=="Low pleiotropy")$gene)
high_pleiotropy_genes = unique(filter(eQTL_dataset, Pleiotropy_category=="High pleiotropy")$gene)

#----------------------------#
# 4) Run enrichment analyses #
#----------------------------#

### Parallelism of gene pools ###
pval1 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", significant_sg_pool, 100000)
pval2 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", highest_sg_pool, 100000)
pval3 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", evol_response_pool, 100000)
pval4 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", common_pool, 100000)
pval5 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", hub_genes, 100000)

### Parallelism of eQTLs ###
pval6  = run_gene_permutation_tests(gene_dataset, "isFullyParallel", eQTL_carriers, 100000)
pval7  = run_gene_permutation_tests(gene_dataset, "isFullyParallel", eQTL_phenotypes, 100000)
pval8  = run_gene_permutation_tests(gene_dataset, "isFullyParallel", no_pleiotropy_genes, 100000)
pval9  = run_gene_permutation_tests(gene_dataset, "isFullyParallel", low_pleiotropy_genes, 100000)
pval10 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", high_pleiotropy_genes, 100000)

### Cross comparisons ###
pval11 = run_gene_permutation_tests(gene_dataset, "isFullyParallel", hub_genes, 100000)
pval12 = run_gene_permutation_tests(gene_dataset, "significant_sg_pool", hub_genes, 100000)
pval13 = run_gene_permutation_tests(gene_dataset, "eQTL_carrier", hub_genes, 100000)
pval14 = run_gene_permutation_tests(gene_dataset, "eQTL_phenotype", hub_genes, 100000)

#-------------------------------#
# 5) Save the results           #
#-------------------------------#
D = data.frame(
  c("significant_sg_pool", "highest_sg_pool", "evol_response_pool", "common_pool", "hub_genes", "eQTL_carriers", "eQTL_phenotypes", "no_pleiotropy_genes", "low_pleiotropy_genes", "high_pleiotropy_genes"),
  c(pval1, pval2, pval3, pval4, pval5, pval6, pval7, pval8, pval9, pval10)
)
names(D) = c("gene_pool", "pvalue")
write.table(D, "./analysis/direct_indirect_selection_analysis/data/gene_parallelism_permutations.txt", row.names=F, col.names=T, quote=F)

D = data.frame(
  c("hub_genes", "hub_genes", "hub_genes", "hub_genes"),
  c("isFullyParallel", "significant_sg_pool", "eQTL_carrier", "eQTL_phenotype"),
  c(pval11, pval12, pval13, pval14)
)
names(D) = c("gene_pool", "enriched_in", "pvalue")
write.table(D, "./analysis/direct_indirect_selection_analysis/data/gene_cross_permutations.txt", row.names=F, col.names=T, quote=F)


