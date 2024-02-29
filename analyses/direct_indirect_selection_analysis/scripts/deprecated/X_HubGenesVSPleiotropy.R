#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 6_HubGenesVSPleiotropy.R
# ------------------------
# Compare eQTL pleiotropy with WGCNA module size.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("cowplot")
library("ggpubr")
library("RColorBrewer")


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#-----------------------------#
# 1) Load the data            #
#-----------------------------#
SNP_dataset            = readRDS("./analyses/direct_indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset           = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_dataset.rds")
eQTL_phenotype_dataset = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_phenotype_dataset.rds")
gene_connectivity      = read.table("./data/tribolium_modules/results_eva/Connectivity_HD.txt", sep="\t", h=T)
rownames(gene_dataset) = gene_dataset$gene

gene_connectivity          = gene_connectivity[order(gene_connectivity$Module, gene_connectivity$kWithin, decreasing=T),]
gene_dataset$module_number = rep(0, dim(gene_dataset)[1])
module_counter             = 1
for(module in unique(gene_connectivity$Module))
{
  N                                       = length(rownames(filter(gene_connectivity, Module==module)))
  gene_list                               = unique(rownames(filter(gene_connectivity, Module==module))[1:N])
  gene_dataset[gene_dataset$gene%in%gene_list,"module_number"] = module_counter
  module_counter                          = module_counter+1
}

D = data.frame()
for (module in unique(gene_dataset$module_number))
{
  hub_genes = filter(gene_dataset, module_number==module & hub_gene==1)$gene
  print(mean(filter(eQTL_dataset, gene%in%hub_genes)$Pleiotropy))
}

