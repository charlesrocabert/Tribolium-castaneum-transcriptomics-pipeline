#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_Statistics.R
# --------------
# Collect some statistics for the manuscript.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")

### Return the stars corresponding to pvalues ###
get_pvalue_stars <- function( pvalue )
{
  star = "ns"
  if (pvalue < 0.0001)
  {
    star = "****"
  }
  if (pvalue >= 0.0001 & pvalue < 0.001)
  {
    star = "***"
  }
  if (pvalue >= 0.001 & pvalue < 0.01)
  {
    star = "**"
  }
  if (pvalue >= 0.01 & pvalue < 0.05)
  {
    star = "*"
  }
  if (pvalue >= 0.05 & pvalue < 0.1)
  {
    star = "."
  }
  if (pvalue >= 0.1)
  {
    star = "ns"
  }
  return(star)
}

### Aggregator function ###
aggregator_function <- function( x, aggregator )
{
  if (aggregator == "sum")
  {
    return(sum(x, na.rm=T))
  }
  if (aggregator == "mean")
  {
    return(mean(x, na.rm=T))
  }
  if (aggregator == "median")
  {
    return(median(x, na.rm=T))
  }
}

### Run one permutations test ###
run_permutations_test <- function( dataset, target_variable, source_variable, source_list, nb_reps, aggregator )
{
  x              = dataset[,target_variable]
  pos            = which(dataset[,source_variable]%in%source_list)
  original_score = aggregator_function(x[pos], aggregator)
  permutations   = t(sapply(1:nb_reps,function(i){
    score = aggregator_function(sample(x)[pos], aggregator)
    return(score)
  }))
  qval = sum(permutations<original_score)/nb_reps
  return(1-qval)
}

### Run a grid of permutation tests ###
run_grid_permutations_test <- function( dataset, target_variables, source_classes, source_variable, nb_reps, aggregator )
{
  N   = length(target_variables)
  M   = length(source_classes)
  RES = matrix(rep(1.0, N*M), ncol=M)
  for(i in seq(1,N))
  {
    for(j in seq(1,M))
    {
      target_var   = target_variables[i]
      source_class = source_classes[j]
      source_list  = dataset[dataset[,source_class]==1,source_variable]
      pval = run_permutations_test(dataset, target_var, source_variable, source_list, nb_reps, aggregator)
      RES[i,j] = pval
      print(paste(target_var, source_class, pval))
    }
  }
  rownames(RES) = target_variables
  colnames(RES) = source_classes
  return(RES)
}

##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#------------------#
# 1) Load the data #
#------------------#
SNP_dataset            = readRDS("./analyses/direct_indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_dataset.rds")
eQTL_phenotype_dataset = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")
eQTL_carrier_dataset   = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_carrier_dataset.rds")

#----------------------#
# 2) Permutation tests #
#----------------------#
target_variables = c("isShifting", "isPartiallyParallel", "isFullyParallel")
source_classes   = c("significant_sg_pool", "highest_sg_pool", "significant_de_pool", "highest_de_pool")
grid_res         = run_grid_permutations_test(SNP_dataset, target_variables, source_classes, "ID", 10000, "sum")
grid_res
write.table(grid_res, "./analyses/direct_indirect_selection_analysis/data/grid_permutation_tests.csv", row.names=T, col.names=T, sep="\t", quote=F)
stop()
grid_stars = grid_res
for(i in seq(1,dim(grid_res)[1]))
{
  for(j in seq(1,dim(grid_res)[2]))
  {
    grid_stars[i,j] = get_pvalue_stars(grid_res[i,j])
  }
}
grid_res
dim(filter(SNP_dataset, Parallelism>1))
mean(filter(SNP_dataset, significant_sg_pool==1)$Parallelism)
mean(SNP_dataset$Parallelism)

filter(SNP_dataset, highest_de_pool==1)
gene_list = filter(gene_dataset, eQTL_carrier==1)$gene
run_permutations_test(gene_dataset, "isFullyParallel", "gene", gene_list, 10000, "sum")

#----------------------#
# 3) General counts    #
#----------------------#
dim(SNP_dataset)
length(unique(SNP_dataset$gene))
mean(table(SNP_dataset$gene))

length(unique(eQTL_dataset$ID))
length(unique(eQTL_dataset$gene))

sum(SNP_dataset$Parallelism==1)
sum(SNP_dataset$Parallelism>=5)
sum(SNP_dataset$Parallelism>1 & SNP_dataset$Parallelism<5)

sum(gene_dataset$Parallelism==1)
sum(gene_dataset$Parallelism>=5)
sum(gene_dataset$Parallelism>1 & gene_dataset$Parallelism<5)

sum(SNP_dataset$significant_sg_pool)
sum(gene_dataset$significant_sg_pool)

sum(gene_dataset$eQTL_carrier)
sum(gene_dataset$eQTL_phenotype)

sum(gene_dataset$hub_gene)

range(gene_dataset$Pleiotropy)

#----------------------#
# 4) Other tests       #
#----------------------#
gene_list = filter(gene_dataset, eQTL_carrier==1)$gene
run_permutations_test(gene_dataset, "isFullyParallel", "gene", gene_list, 10000, "sum")

gene_list = filter(gene_dataset, eQTL_phenotype==1)$gene
run_permutations_test(gene_dataset, "isFullyParallel", "gene", gene_list, 10000, "sum")

gene_list = filter(gene_dataset, hub_gene==1)$gene
run_permutations_test(gene_dataset, "isFullyParallel", "gene", gene_list, 10000, "sum")

gene_list = filter(gene_dataset, hub_gene==1)$gene
run_permutations_test(gene_dataset, "eQTL_carrier", "gene", gene_list, 10000, "sum")

fisher.test(gene_dataset$hub_gene, gene_dataset$significant_sg_pool)

fisher.test(gene_dataset$hub_gene, gene_dataset$eQTL_carrier)
fisher.test(SNP_dataset$highest_de_pool, SNP_dataset$is_ASE)

x1 = filter(gene_dataset, is_ASE==0)$abs_logFC
x2 = filter(gene_dataset, is_ASE==1)$abs_logFC
wilcox.test(x1,x2)
aggregator_function(x1, "mean")
aggregator_function(x2, "mean")


x1 = filter(gene_dataset, eQTL_phenotype==0)$abs_logFC
x2 = filter(gene_dataset, eQTL_phenotype==1)$abs_logFC
wilcox.test(x1,x2)
aggregator_function(x1, "mean")
aggregator_function(x2, "mean")

x1 = filter(gene_dataset, hub_gene==0)$abs_logFC
x2 = filter(gene_dataset, hub_gene==1)$abs_logFC
wilcox.test(x1,x2)
aggregator_function(x1, "mean")
aggregator_function(x2, "mean")

x1 = filter(gene_dataset, eQTL_carrier==1)$abs_logFC
x2 = filter(gene_dataset, eQTL_phenotype==1)$abs_logFC
wilcox.test(x1,x2)
aggregator_function(x1, "mean")
aggregator_function(x2, "mean")


gene_list = filter(gene_dataset, Parallelism>0)$gene
run_permutations_test(gene_dataset, "abs_logFC", "gene", gene_list, 10000, "mean")
head(gene_dataset)
mean(x1, na.rm=T)
mean(x2, na.rm=T)
plot(gene_dataset$abs_logFC, gene_dataset$log_net_selection)
