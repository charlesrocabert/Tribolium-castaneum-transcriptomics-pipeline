#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_Statistics.R
# --------------
# Show some statistics on the datasets.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("cowplot")
library("rbioapi")

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

### Plot parallelism distribution as barplots ###
plot_parallel_distribution <- function( eQTL_dataset, SNP_dataset, gene_dataset )
{
  h1 = ggplot(SNP_dataset, aes(Parallelism)) +
    geom_bar() +
    scale_x_continuous(breaks=0:7) +
    xlab("Parallelism") +
    ylab("Count") +
    ggtitle("SNPs") +
    theme_classic() +
    geom_text(aes(label=..count.., y=..count..), stat="count", vjust=-.5, size=3)
  h2 = ggplot(eQTL_dataset, aes(Parallelism)) +
    geom_bar() +
    scale_x_continuous(breaks=0:7) +
    xlab("Parallelism") +
    ylab("Count") +
    ggtitle("eQTLs") +
    theme_classic() +
    geom_text(aes(label=..count.., y=..count..), stat="count", vjust=-.5, size=3)
  h3 = ggplot(gene_dataset, aes(Parallelism)) +
    geom_bar() +
    scale_x_continuous(breaks=0:7) +
    xlab("Parallelism") +
    ylab("Count") +
    ggtitle("Genes") +
    theme_classic() +
    geom_text(aes(label=..count.., y=..count..), stat="count", vjust=-.5, size=3)
  p = plot_grid(h1, h3, ncol=2, labels="AUTO")
  return(p)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

eQTLs        = load_eQTLs("EXPRESSION")
SNP_dataset  = readRDS("./analyses/direct_indirect_selection_analysis/data/SNP_dataset.rds")
gene_dataset = readRDS("./analyses/direct_indirect_selection_analysis/data/gene_dataset.rds")
eQTL_dataset = readRDS("./analyses/direct_indirect_selection_analysis/data/eQTL_dataset.rds")

#----------------------------------------------#
# 0) Save significant shifts data              #
#----------------------------------------------#
to_keep = c("gene", "ID", "POS", "CHROM", "L1_shift", "L2_shift", "L3_shift", "L5_shift", "L6_shift", "Mx1_shift", "Mx2_shift", "Parallelism", "is_ASE")
write.table(SNP_dataset[,to_keep], file="./analyses/direct_indirect_selection_analysis/manuscript/DataS8/DataS8.csv", sep=";", col.names=T, row.names=F, quote=F)

#----------------------------------------------#
# 1) Get some counts                           #
#----------------------------------------------#

### All SNPs ###
n = sum(SNP_dataset$Parallelism==1)
p = n/dim(SNP_dataset)[1]*100
print(paste0("> ",n," single shifting SNPs (",round(p,2),"%)"))

n = sum(SNP_dataset$Parallelism%in%c(2,3,4))
p = n/dim(SNP_dataset)[1]*100
print(paste0("> ",n," partially parallel SNPs (",round(p,2),"%)"))

n = sum(SNP_dataset$Parallelism>4)
p = n/dim(SNP_dataset)[1]*100
print(paste0("> ",n," fully parallel SNPs (",round(p,2),"%)"))

### Genes ###
n = sum(gene_dataset$Parallelism==1)
p = n/dim(gene_dataset)[1]*100
print(paste0("> ",n," single shifting genes (",round(p,2),"%)"))

n = sum(gene_dataset$Parallelism%in%c(2,3,4))
p = n/dim(gene_dataset)[1]*100
print(paste0("> ",n," partially parallel genes (",round(p,2),"%)"))

n = sum(gene_dataset$Parallelism>4)
p = n/dim(gene_dataset)[1]*100
print(paste0("> ",n," fully parallel genes (",round(p,2),"%)"))

### eQTLs ###
n = sum(eQTL_dataset$Parallelism==1)
p = n/dim(eQTL_dataset)[1]*100
print(paste0("> ",n," single shifting eQTLs (",round(p,2),"%)"))

n = sum(eQTL_dataset$Parallelism%in%c(2,3,4))
p = n/dim(eQTL_dataset)[1]*100
print(paste0("> ",n," partially parallel eQTLs (",round(p,2),"%)"))

n = sum(eQTL_dataset$Parallelism>4)
p = n/dim(eQTL_dataset)[1]*100
print(paste0("> ",n," fully parallel eQTLs (",round(p,2),"%)"))

#----------------------------------------------#
# 2) Parallelism barplots                      #
#----------------------------------------------#
p = plot_parallel_distribution(eQTL_dataset, SNP_dataset, gene_dataset)
ggsave("analysis/direct_indirect_selection_analysis/plots/parallelism_barplots.pdf", p, width=8, height=4, units="in")

#----------------------------------------------#
# 3) Save various gene lists for GO enrichment #
#----------------------------------------------#
significant_sg_genes = filter(gene_dataset, significant_sg_pool==1)$gene
parallel_genes       = filter(gene_dataset, isFullyParallel==1)$gene
fully_parallel_genes = filter(gene_dataset, Parallelism==7)$gene
eQTL_carriers        = filter(gene_dataset, eQTL_carrier==1)$gene
eQTL_phenotypes      = filter(gene_dataset, eQTL_phenotype==1)$gene
hub_genes            = filter(gene_dataset, hub_gene==1)$gene

write.table(significant_sg_genes, "analysis/direct_indirect_selection_analysis/data/significant_sg_genes.txt", quote=F, row.names=F, col.names=F)
write.table(parallel_genes, "analysis/direct_indirect_selection_analysis/data/parallel_genes.txt", quote=F, row.names=F, col.names=F)
write.table(fully_parallel_genes, "analysis/direct_indirect_selection_analysis/data/fully_parallel_genes.txt", quote=F, row.names=F, col.names=F)
write.table(eQTL_carriers, "analysis/direct_indirect_selection_analysis/data/eQTL_carriers.txt", quote=F, row.names=F, col.names=F)
write.table(eQTL_phenotypes, "analysis/direct_indirect_selection_analysis/data/eQTL_phenotypes.txt", quote=F, row.names=F, col.names=F)
write.table(hub_genes, "analysis/direct_indirect_selection_analysis/data/hub_genes.txt", quote=F, row.names=F, col.names=F)

#----------------------------------------------#
# 4) Sg genes statistics                       #
#----------------------------------------------#
n = sum(gene_dataset$significant_sg_pool==1)
p = n/dim(gene_dataset)[1]*100
print(paste0("> ",n," significant Sg genes (",round(p,2),"%)"))

n = sum(SNP_dataset$significant_sg_pool==1)
p = n/dim(SNP_dataset)[1]*100
print(paste0("> ",n," significant Sg SNPs (",round(p,2),"%)"))

n = dim(filter(eQTLs, phenotype%in%filter(gene_dataset, significant_sg_pool==1)$gene))[1]
p = n/dim(eQTLs)[1]*100
print(paste0("> ",n," significant Sg eQTLs (",round(p,2),"%)"))

