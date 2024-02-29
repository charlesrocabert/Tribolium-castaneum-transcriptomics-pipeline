#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_BuildDatasets.R
# -----------------
# Build the genetic map dataset prior to analysis.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")

### Mean wrapper function ###
my_mean <- function( x )
{
  return(mean(x, na.rm=T))
}

### Sum wrapper function ###
my_sum <- function( x )
{
  return(sum(x, na.rm=T))
}

### Plot SNP dataset stats ###
SNP_dataset_statistics <- function( SNP_dataset )
{
  nb_snps      = dim(SNP_dataset)[1]
  nb_genes     = length(unique(SNP_dataset$gene))
  nb_sg_snps   = sum(SNP_dataset$significant_sg_pool==1)
  nb_hub_snps  = sum(SNP_dataset$hub_gene==1)
  nb_sg_genes  = sum(SNP_dataset[!duplicated(SNP_dataset$gene),]$significant_sg_pool==1)
  nb_hub_genes = sum(SNP_dataset[!duplicated(SNP_dataset$gene),]$hub_gene==1)
  l1 = paste0("> nb genes = ", nb_genes, " (", nb_snps, " SNPs)")
  l2 = paste0("> nb Sg genes = ", nb_sg_genes, " (", nb_sg_snps, " SNPs)")
  l3 = paste0("> nb hub genes = ", nb_hub_genes, " (", nb_hub_snps, " SNPs)")
  cat(paste0(l1, "\n", l2, "\n", l3, "\n"))
}

### Plot SNP dataset stats ###
gene_dataset_statistics <- function( gene_dataset )
{
  nb_genes     = length(gene_dataset$gene)
  nb_sg_genes  = sum(gene_dataset$significant_sg_pool==1)
  nb_hub_genes = sum(gene_dataset$hub_gene==1)
  l1 = paste0("> nb genes = ", nb_genes)
  l2 = paste0("> nb Sg genes = ", nb_sg_genes)
  l3 = paste0("> nb hub genes = ", nb_hub_genes)
  cat(paste0(l1, "\n", l2, "\n", l3, "\n"))
}

### Pivot HD lines AFC values ###
pivot_AFC_HD_table <- function( AFC_HD, transform_AFCs )
{
  #------------------------------------------#
  # 1) Load L1 data                          #
  #------------------------------------------#
  L1_G1_HD     = filter(AFC_HD, LINE == "L1")[,c("ID","AF_G1")]
  L1_G21_HD    = filter(AFC_HD, LINE == "L1")[,c("ID","AF_G21")]
  L1_HD        = merge(L1_G1_HD, L1_G21_HD, by="ID")
  L1_HD$L1_AFC = L1_HD$AF_G21-L1_HD$AF_G1
  names(L1_HD) = c("ID", "L1_HD_G1", "L1_HD_G21", "L1_AFC")
  #------------------------------------------#
  # 2) Load L2 data                          #
  #------------------------------------------#
  L2_G1_HD     = filter(AFC_HD, LINE == "L2")[,c("ID","AF_G1")]
  L2_G21_HD    = filter(AFC_HD, LINE == "L2")[,c("ID","AF_G21")]
  L2_HD        = merge(L2_G1_HD, L2_G21_HD, by="ID")
  L2_HD$L2_AFC = L2_HD$AF_G21-L2_HD$AF_G1
  names(L2_HD) = c("ID", "L2_HD_G1", "L2_HD_G21", "L2_AFC")
  #------------------------------------------#
  # 3) Load L3 data                          #
  #------------------------------------------#
  L3_G1_HD     = filter(AFC_HD, LINE == "L3")[,c("ID","AF_G1")]
  L3_G21_HD    = filter(AFC_HD, LINE == "L3")[,c("ID","AF_G21")]
  L3_HD        = merge(L3_G1_HD, L3_G21_HD, by="ID")
  L3_HD$L3_AFC = L3_HD$AF_G21-L3_HD$AF_G1
  names(L3_HD) = c("ID", "L3_HD_G1", "L3_HD_G21", "L3_AFC")
  #------------------------------------------#
  # 4) Load L4 data                          #
  #------------------------------------------#
  L5_G1_HD     = filter(AFC_HD, LINE == "L5")[,c("ID","AF_G1")]
  L5_G21_HD    = filter(AFC_HD, LINE == "L5")[,c("ID","AF_G21")]
  L5_HD        = merge(L5_G1_HD, L5_G21_HD, by="ID")
  L5_HD$L5_AFC = L5_HD$AF_G21-L5_HD$AF_G1
  names(L5_HD) = c("ID", "L5_HD_G1", "L5_HD_G21", "L5_AFC")
  #------------------------------------------#
  # 5) Load L5 data                          #
  #------------------------------------------#
  L6_G1_HD     = filter(AFC_HD, LINE == "L6")[,c("ID","AF_G1")]
  L6_G21_HD    = filter(AFC_HD, LINE == "L6")[,c("ID","AF_G21")]
  L6_HD        = merge(L6_G1_HD, L6_G21_HD, by="ID")
  L6_HD$L6_AFC = L6_HD$AF_G21-L6_HD$AF_G1
  names(L6_HD) = c("ID", "L6_HD_G1", "L6_HD_G21", "L6_AFC")
  #-------------------------------------------#
  # 6) Load Mx1 data                          #
  #-------------------------------------------#
  Mx1_G1_HD      = filter(AFC_HD, LINE == "Mx1")[,c("ID","AF_G1")]
  Mx1_G21_HD     = filter(AFC_HD, LINE == "Mx1")[,c("ID","AF_G21")]
  Mx1_HD         = merge(Mx1_G1_HD, Mx1_G21_HD, by="ID")
  Mx1_HD$Mx1_AFC = Mx1_HD$AF_G21-Mx1_HD$AF_G1
  names(Mx1_HD) = c("ID", "Mx1_HD_G1", "Mx1_HD_G21", "Mx1_AFC")
  #-------------------------------------------#
  # 7) Load Mx2 data                          #
  #-------------------------------------------#
  Mx2_G1_HD      = filter(AFC_HD, LINE == "Mx2")[,c("ID","AF_G1")]
  Mx2_G21_HD     = filter(AFC_HD, LINE == "Mx2")[,c("ID","AF_G21")]
  Mx2_HD         = merge(Mx2_G1_HD, Mx2_G21_HD, by="ID")
  Mx2_HD$Mx2_AFC = Mx2_HD$AF_G21-Mx2_HD$AF_G1
  names(Mx2_HD) = c("ID", "Mx2_HD_G1", "Mx2_HD_G21", "Mx2_AFC")
  #-------------------------------------------#
  # 8) Build the AFC table                    #
  #-------------------------------------------#
  AFC_TABLE = merge(L1_HD, L2_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L3_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L5_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, L6_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx1_HD, by="ID")
  AFC_TABLE = merge(AFC_TABLE, Mx2_HD, by="ID")
  #-------------------------------------------#
  # 9) Normalize AFs                          #
  #-------------------------------------------#
  coln = c("L1_HD_G1", "L1_HD_G21", "L2_HD_G1", "L2_HD_G21", "L3_HD_G1", "L3_HD_G21", "L5_HD_G1", "L5_HD_G21", "L6_HD_G1", "L6_HD_G21", "Mx1_HD_G1", "Mx1_HD_G21", "Mx2_HD_G1", "Mx2_HD_G21")
  if (transform_AFCs)
  {
    AFC_TABLE[,coln] = asin(sqrt(AFC_TABLE[,coln]))
    AFC_TABLE[,coln] = scale(AFC_TABLE[,coln], center=T, scale=T)
  }
  #-------------------------------------------#
  # 10) Split ID in pos and chr               #
  #-------------------------------------------#
  out  = strsplit(AFC_TABLE$ID,'-') 
  out2 = do.call(rbind, out)
  AFC_TABLE$CHR = out2[,1]
  AFC_TABLE$POS = out2[,2]
  #-------------------------------------------#
  # 11) Save absolute values of AFCs per line #
  #-------------------------------------------#
  AFC_TABLE$L1_AFC_abs  = abs(AFC_TABLE$L1_AFC)
  AFC_TABLE$L2_AFC_abs  = abs(AFC_TABLE$L2_AFC)
  AFC_TABLE$L3_AFC_abs  = abs(AFC_TABLE$L3_AFC)
  AFC_TABLE$L5_AFC_abs  = abs(AFC_TABLE$L5_AFC)
  AFC_TABLE$L6_AFC_abs  = abs(AFC_TABLE$L6_AFC)
  AFC_TABLE$Mx1_AFC_abs = abs(AFC_TABLE$Mx1_AFC)
  AFC_TABLE$Mx2_AFC_abs = abs(AFC_TABLE$Mx2_AFC)
  #-------------------------------------------#
  # 12) Return the table                      #
  #-------------------------------------------#
  return(AFC_TABLE)
}

### Load eQTLs ###
load_eQTLs <- function( phenotype )
{
  #--------------------------#
  # 1) Load eQTLs data       #
  #--------------------------#
  EQTLs    = readRDS(paste0("./data/tribolium_eqtl/significant/HD_G1_Tcas3.30_imputed_",phenotype,"_significant.rds"))
  EQTLs$ID = EQTLs$rs
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

### Build the SNP dataset ###
build_SNP_dataset <- function()
{
  #-------------------#
  # 1) Load datasets  #
  #-------------------#
  
  ### Load allelic frequency changes (AFCs) data ###
  AFC_HD = readRDS(paste0("./data/tribolium_afc/AFC_HD.rds"))
  AFC_HD = pivot_AFC_HD_table(AFC_HD, FALSE)
  AFC_HD = AFC_HD[,c("ID", "POS", "L1_AFC", "L2_AFC", "L3_AFC", "L5_AFC", "L6_AFC", "Mx1_AFC", "Mx2_AFC", "L1_AFC_abs", "L2_AFC_abs", "L3_AFC_abs", "L5_AFC_abs", "L6_AFC_abs", "Mx1_AFC_abs", "Mx2_AFC_abs")]
  AFC_HD = AFC_HD[!duplicated(AFC_HD$ID),]
  
  ### Load SNP parallelism information ###
  L1_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L1.rds")
  L2_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L2.rds")
  L3_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L3.rds")
  L5_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L5.rds")
  L6_SIGNIF  = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_L6.rds")
  Mx1_SIGNIF = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_Mx1.rds")
  Mx2_SIGNIF = readRDS("./data/tribolium_afc/signifSNP_AFC/signifSNP_AFC_Mx2.rds")
  
  L1_SIGNIF  = as.data.frame(L1_SIGNIF$alpha05)
  L2_SIGNIF  = as.data.frame(L2_SIGNIF$alpha05)
  L3_SIGNIF  = as.data.frame(L3_SIGNIF$alpha05)
  L5_SIGNIF  = as.data.frame(L5_SIGNIF$alpha05)
  L6_SIGNIF  = as.data.frame(L6_SIGNIF$alpha05)
  Mx1_SIGNIF = as.data.frame(Mx1_SIGNIF$alpha05)
  Mx2_SIGNIF = as.data.frame(Mx2_SIGNIF$alpha05)
  
  names(L1_SIGNIF)  = c("index", "pvalue")
  names(L2_SIGNIF)  = c("index", "pvalue")
  names(L3_SIGNIF)  = c("index", "pvalue")
  names(L5_SIGNIF)  = c("index", "pvalue")
  names(L6_SIGNIF)  = c("index", "pvalue")
  names(Mx1_SIGNIF) = c("index", "pvalue")
  names(Mx2_SIGNIF) = c("index", "pvalue")
  
  L1_SIGNIF$index   = as.numeric(L1_SIGNIF$index)
  L1_SIGNIF$pvalue  = as.numeric(L1_SIGNIF$pvalue)
  L2_SIGNIF$index   = as.numeric(L2_SIGNIF$index)
  L2_SIGNIF$pvalue  = as.numeric(L2_SIGNIF$pvalue)
  L3_SIGNIF$index   = as.numeric(L3_SIGNIF$index)
  L3_SIGNIF$pvalue  = as.numeric(L3_SIGNIF$pvalue)
  L5_SIGNIF$index   = as.numeric(L5_SIGNIF$index)
  L5_SIGNIF$pvalue  = as.numeric(L5_SIGNIF$pvalue)
  L6_SIGNIF$index   = as.numeric(L6_SIGNIF$index)
  L6_SIGNIF$pvalue  = as.numeric(L6_SIGNIF$pvalue)
  Mx1_SIGNIF$index  = as.numeric(Mx1_SIGNIF$index)
  Mx1_SIGNIF$pvalue = as.numeric(Mx1_SIGNIF$pvalue)
  Mx2_SIGNIF$index  = as.numeric(Mx2_SIGNIF$index)
  Mx2_SIGNIF$pvalue = as.numeric(Mx2_SIGNIF$pvalue)
  
  ### Load SNP annotation ###
  ANNOTATION        = read.table("./data/tribolium_snp/snp_table_ALL_Tcas3.30_raw_SNP.csv", h=T, sep="\t")
  ANNOTATION$gene   = ANNOTATION$Feature_id
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "Annotation", "Putative_impact", "gene"))
  
  ### Load list of genes under selection ###
  
  # 1) Genes with a significant s_g AND Genes with the highest abs(s_g)
  # -------------------------------------------------------------------
  GENE_FITNESS_COR_DATASET          = read.table("./data/experiment_data/GeneticCor_Fitness_HD.csv", sep=";", h=T, dec=",")
  GENE_FITNESS_COR_DATASET$Mean_abs = abs(GENE_FITNESS_COR_DATASET$Mean)
  GENE_FITNESS_COR_DATASET          = GENE_FITNESS_COR_DATASET[order(GENE_FITNESS_COR_DATASET$Mean_abs, decreasing=T),]
  significant_sg_pool               = filter(GENE_FITNESS_COR_DATASET, sig>0)$gene
  highest_sg_pool                   = GENE_FITNESS_COR_DATASET$gene[1:106]
  
  # 2) Genes showing an evolutionary response (Koch & Guillaume, 2020)
  # ------------------------------------------------------------------
  EVOL_RESPONSE_DATASET           = read.table("./data/experiment_data/Evol-Diff-Expr_HD_G22.txt", sep="\t", h=T)
  EVOL_RESPONSE_DATASET$abs_logFC = abs(EVOL_RESPONSE_DATASET$logFC)
  EVOL_RESPONSE_DATASET           = EVOL_RESPONSE_DATASET[order(EVOL_RESPONSE_DATASET$abs_logFC, decreasing=T),]
  significant_de_pool             = rownames(filter(EVOL_RESPONSE_DATASET, DE!=0))
  highest_de_pool                 = rownames(EVOL_RESPONSE_DATASET[1:135,])
  
  # 3) Create a common group for all previous sets
  # ----------------------------------------------
  common_pool = c(significant_sg_pool, highest_sg_pool, significant_de_pool, highest_de_pool)
  
  # 4) Genes in expression modules correlated to fitness
  # ----------------------------------------------------
  EXPRESSION_MODULES      = read.table("./data/tribolium_modules/HD_G1_Tcas3.30_EXPRESSION_significant_modules.txt", sep="\t", h=T)
  expression_modules_pool = EXPRESSION_MODULES$gene
  
  # 5) WGCNA hub genes
  # ------------------
  gene_connectivity = read.table("./data/tribolium_modules/results_eva/Connectivity_HD.txt", sep="\t", h=T)
  gene_connectivity = gene_connectivity[order(gene_connectivity$Module, gene_connectivity$kWithin, decreasing=T),]
  hub_genes         = c()
  module_genes      = c()
  for(module in unique(gene_connectivity$Module))
  {
    N            = length(rownames(filter(gene_connectivity, Module==module)))
    hub_genes    = c(hub_genes, rownames(filter(gene_connectivity, Module==module))[1:10])
    module_genes = c(module_genes, rownames(filter(gene_connectivity, Module==module))[11:N])
  }
  
  ### Load phenotype positions ###
  PHENOTYPE_POSITIONS = read.table("./data/tribolium_eqtl/gene_pos_Tcas3.30.csv", sep="\t", h=T)
  names(PHENOTYPE_POSITIONS) = c("gene", "gene_chr", "gene_start", "gene_end")
  
  ### Load ASE ###
  #  "Several additional filtering conditions for ASE SNPs were set at the
  #  population level: at least 30 heterozygotes and at least 10 heterozygotes
  #  that displayed ASE" (Liu et al. 2020)
  HET_SNPS         = read.table("./data/tribolium_ase/HD_G1_Tcas3.30_imputed_HET_SNPs.txt", h=T, row.names=1, check.names=F)
  ASE_SNPS         = read.table("./data/tribolium_ase/HD_G1_Tcas3.30_imputed_SIGNIF_SNPs.txt", h=T, row.names=1, check.names=F)
  HET_COUNT        = data.frame(names(rowSums(HET_SNPS)), rowSums(HET_SNPS))
  names(HET_COUNT) = c("ID", "HET_COUNT")
  ASE_COUNT        = data.frame(names(rowSums(ASE_SNPS)), rowSums(ASE_SNPS))
  names(ASE_COUNT) = c("ID", "ASE_COUNT")
  ASE_DATA         = merge(HET_COUNT, ASE_COUNT, by="ID")
  ASE_DATA$IS_ASE  = as.numeric((ASE_DATA$HET_COUNT >= 30 & ASE_DATA$ASE_COUNT >= 10))
  
  ### Load fitness covariance results ###
  FITNESS_DATA           = read.table("./data/experiment_data/Selection-Cov-Beta-HD.csv", h=T, sep=",")
  rownames(FITNESS_DATA) = FITNESS_DATA$gene
  
  #-------------------#
  # 2) Merge datasets #
  #-------------------#
  
  ### Create the final dataset ###
  DATA = data.frame(AFC_HD)
  
  ### Merge parallelism information ###
  DATA$L1_shift                   = as.numeric(DATA$ID%in%DATA$ID[L1_SIGNIF$index])
  DATA$L2_shift                   = as.numeric(DATA$ID%in%DATA$ID[L2_SIGNIF$index])
  DATA$L3_shift                   = as.numeric(DATA$ID%in%DATA$ID[L3_SIGNIF$index])
  DATA$L5_shift                   = as.numeric(DATA$ID%in%DATA$ID[L5_SIGNIF$index])
  DATA$L6_shift                   = as.numeric(DATA$ID%in%DATA$ID[L6_SIGNIF$index])
  DATA$Mx1_shift                  = as.numeric(DATA$ID%in%DATA$ID[Mx1_SIGNIF$index])
  DATA$Mx2_shift                  = as.numeric(DATA$ID%in%DATA$ID[Mx2_SIGNIF$index])
  DATA$Parallelism                = DATA$L1_shift+DATA$L2_shift+DATA$L3_shift+DATA$L5_shift+DATA$L6_shift+DATA$Mx1_shift+DATA$Mx2_shift
  DATA$isShifting                 = as.numeric(DATA$Parallelism>0)
  DATA$isSingleShifting           = as.numeric(DATA$Parallelism==1)
  DATA$isPartiallyorFullyParallel = as.numeric(DATA$Parallelism>1)
  DATA$isPartiallyParallel        = as.numeric(DATA$Parallelism>1 & DATA$Parallelism<5)
  DATA$isFullyParallel            = as.numeric(DATA$Parallelism>4)
  
  ### Build parallel categories ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No significant shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial or full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No significant shift","Single shift","Partial or full parallelism"))
  
  ### Merge annotation ###
  DATA = merge(DATA, ANNOTATION_SELECT, by="ID")
  
  ### Merge with evolutionary response ###
  DATA$logFC     = EVOL_RESPONSE_DATASET[DATA$gene, "logFC"]
  DATA$abs_logFC = EVOL_RESPONSE_DATASET[DATA$gene, "abs_logFC"]
  
  ### Merge with genes under selection ###
  DATA$significant_sg_pool = as.numeric(DATA$gene%in%significant_sg_pool)
  DATA$highest_sg_pool     = as.numeric(DATA$gene%in%highest_sg_pool)
  DATA$significant_de_pool = as.numeric(DATA$gene%in%significant_de_pool)
  DATA$highest_de_pool     = as.numeric(DATA$gene%in%highest_de_pool)
  DATA$common_pool         = as.numeric(DATA$gene%in%common_pool)
  DATA$hub_gene            = as.numeric(DATA$gene%in%hub_genes)
  DATA$module_gene         = as.numeric(DATA$gene%in%module_genes)
  
  ### Merge with gene positions ###
  DATA = merge(DATA, PHENOTYPE_POSITIONS, by="gene")
  
  ### Merge with ASE data ###
  ASE_SNPs    = filter(ASE_DATA, IS_ASE==1)$ID
  DATA$is_ASE = as.numeric(DATA$ID%in%ASE_SNPs)
  
  ### Merge with fitness data ###
  # COV_g  = direct + undirect (net)
  # beta_g = direct only
  #
  DATA$net_selection        = FITNESS_DATA[DATA$gene,"Cov.Mean"]
  DATA$direct_selection     = FITNESS_DATA[DATA$gene,"bgHD.Mean"]
  DATA$diff_selection       = DATA$net_selection-DATA$direct_selection
  DATA$log_net_selection    = log10(abs(DATA$net_selection))
  DATA$log_direct_selection = log10(abs(DATA$direct_selection))
  
  ### Clear memory ###
  rm(AFC_HD)
  rm(L1_SIGNIF)
  rm(L2_SIGNIF)
  rm(L3_SIGNIF)
  rm(L5_SIGNIF)
  rm(L6_SIGNIF)
  rm(Mx1_SIGNIF)
  rm(Mx2_SIGNIF)
  rm(ANNOTATION)
  rm(ANNOTATION_SELECT)
  rm(EVOL_RESPONSE_DATASET)
  rm(GENE_FITNESS_COR_DATASET)
  rm(PHENOTYPE_POSITIONS)
  rm(HET_SNPS)
  rm(ASE_SNPS)
  rm(HET_COUNT)
  rm(ASE_COUNT)
  rm(ASE_DATA)
  rm(FITNESS_DATA)
  rm(ASE_SNPs)
  
  ### Return the dataset ###
  rownames(DATA) = DATA$ID
  return(DATA)
}

### Build the gene dataset ###
build_gene_dataset <- function( SNP_dataset, eQTLs )
{
  ### Collapse data to genes ###
  DATA                     = SNP_dataset[!duplicated(SNP_dataset$gene),]
  DATA$Parallelism         = tapply(SNP_dataset$Parallelism, SNP_dataset$gene, max)[DATA$gene]
  DATA$isShifting          = tapply(SNP_dataset$isShifting, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isShifting          = as.numeric(DATA$isShifting>0)
  DATA$isPartiallyParallel = tapply(SNP_dataset$isPartiallyParallel, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isPartiallyParallel = as.numeric(DATA$isPartiallyParallel>0)
  DATA$isFullyParallel     = tapply(SNP_dataset$isFullyParallel, SNP_dataset$gene, sum)[DATA$gene]
  DATA$isFullyParallel     = as.numeric(DATA$isFullyParallel>0)
  
  ### Add eQTL relationships ###
  DATA$eQTL_carrier          = as.numeric(DATA$gene%in%eQTLs$gene)
  DATA$eQTL_carrier_category = as.character(DATA$eQTL_carrier)
  DATA$eQTL_carrier_category[DATA$eQTL_carrier_category=="0"] = "Other"
  DATA$eQTL_carrier_category[DATA$eQTL_carrier_category=="1"] = "eQTL carrier"
  DATA$eQTL_carrier_category = factor(DATA$eQTL_carrier_category, levels=c("Other","eQTL carrier"))
  
  DATA$eQTL_phenotype          = as.numeric(DATA$gene%in%eQTLs$phenotype)
  DATA$eQTL_phenotype_category = as.character(DATA$eQTL_phenotype)
  DATA$eQTL_phenotype_category[DATA$eQTL_phenotype_category=="0"] = "Other"
  DATA$eQTL_phenotype_category[DATA$eQTL_phenotype_category=="1"] = "eQTL phenotype"
  DATA$eQTL_phenotype_category = factor(DATA$eQTL_phenotype_category, levels=c("Other","eQTL phenotype"))
  
  ### Build parallel categories ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No significant shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial or full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No significant shift","Single shift","Partial or full parallelism"))
  
  ### Return the dataset ###
  rownames(DATA) = DATA$gene
  return(DATA)
}

### Build the eQTL dataset ###
build_eQTL_dataset <- function( SNP_dataset, gene_dataset, eQTLs )
{
  ### Build the initial eQTL dataset ###
  DATA                           = data.frame(eQTLs)
  to_keep                        = names(DATA)[!names(DATA)%in%c("Annotation","Putative_impact","gene")]
  DATA                           = DATA[,to_keep]
  DATA                           = merge(DATA, SNP_dataset, by="ID")
  DATA$net_selection             = gene_dataset[DATA$phenotype, "net_selection"]
  DATA$direct_selection          = gene_dataset[DATA$phenotype, "direct_selection"]
  #DATA                           = drop_na(DATA)
  DATA$weighted_net_selection    = DATA$beta*((DATA$net_selection))
  DATA$weighted_direct_selection = DATA$beta*((DATA$direct_selection))
  
  ### Calculate additional statistics ###
  mean_net_selection             = tapply(DATA$net_selection, DATA$ID, my_mean)
  mean_direct_selection          = tapply(DATA$direct_selection, DATA$ID, my_mean)
  weighted_net_selection         = tapply(DATA$weighted_net_selection, DATA$ID, my_sum)
  weighted_direct_selection      = tapply(DATA$weighted_direct_selection, DATA$ID, my_sum)
  DATA$Pleiotropy                = as.vector(table(DATA$ID)[DATA$ID])
  #DATA                           = DATA[!duplicated(DATA$ID),]
  DATA$mean_net_selection        = mean_net_selection[DATA$ID]
  DATA$mean_direct_selection     = mean_direct_selection[DATA$ID]
  DATA$weighted_net_selection    = weighted_net_selection[DATA$ID]
  DATA$weighted_direct_selection = weighted_direct_selection[DATA$ID]
  DATA$log_net_selection         = log10(abs(DATA$weighted_net_selection))
  DATA$log_direct_selection      = log10(abs(DATA$weighted_direct_selection))
  
  ### Build pleiotropy classes ###
  DATA$Pleiotropy_category = DATA$Pleiotropy
  DATA$Pleiotropy_category[DATA$Pleiotropy_category>1 & DATA$Pleiotropy_category<5] = 2
  DATA$Pleiotropy_category[DATA$Pleiotropy_category>=5] = 3
  DATA$Pleiotropy_category = as.character(DATA$Pleiotropy_category)
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="1"] = "No pleiotropy"
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="2"] = "Low pleiotropy"
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="3"] = "High pleiotropy"
  DATA$Pleiotropy_category = factor(DATA$Pleiotropy_category, levels=c("No pleiotropy","Low pleiotropy","High pleiotropy"))
  
  ### Build parallel classes ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No significant shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial or full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No significant shift","Single shift","Partial or full parallelism"))
  
  ### Clear memory ###
  rm(mean_net_selection)
  rm(mean_direct_selection)
  rm(weighted_net_selection)
  rm(weighted_direct_selection)
  
  ### Return the dataset ###
  #rownames(DATA) = DATA$ID
  return(DATA)
}

### Build the eQTL carrier dataset ###
build_eQTL_carrier_dataset <- function( eQTL_dataset, gene_dataset )
{
  ### Calculate pleiotropy at the gene level ###
  Pleiotropy = c()
  for (g in unique(eQTL_dataset$gene))
  {
    val        = length(unique(filter(eQTL_dataset, gene==g)$phenotype))
    Pleiotropy = c(Pleiotropy, val)
  }
  names(Pleiotropy) = unique(eQTL_dataset$gene)
  
  ### Merge datasets ###
  DATA = data.frame(Pleiotropy, gene_dataset[names(Pleiotropy),"log_net_selection"])
  #DATA = drop_na(DATA)
  
  ### Build pleiotropy classes ###
  # DATA$Pleiotropy_category = DATA$Pleiotropy
  # DATA$Pleiotropy_category[DATA$Pleiotropy_category>1 & DATA$Pleiotropy_category<5] = 2
  # DATA$Pleiotropy_category[DATA$Pleiotropy_category>=5] = 3
  # DATA$Pleiotropy_category = as.character(DATA$Pleiotropy_category)
  # DATA$Pleiotropy_category[DATA$Pleiotropy_category=="1"] = "No pleiotropy"
  # DATA$Pleiotropy_category[DATA$Pleiotropy_category=="2"] = "Low pleiotropy"
  # DATA$Pleiotropy_category[DATA$Pleiotropy_category=="3"] = "High pleiotropy"
  # DATA$Pleiotropy_category = factor(DATA$Pleiotropy_category, levels=c("No pleiotropy","Low pleiotropy", "High pleiotropy"))
  ### Build pleiotropy classes ###
  DATA$Pleiotropy_category = DATA$Pleiotropy
  DATA$Pleiotropy_category[DATA$Pleiotropy_category>1] = 2
  DATA$Pleiotropy_category = as.character(DATA$Pleiotropy_category)
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="1"] = "No pleiotropy"
  DATA$Pleiotropy_category[DATA$Pleiotropy_category=="2"] = "Pleiotropy"
  DATA$Pleiotropy_category = factor(DATA$Pleiotropy_category, levels=c("No pleiotropy","Pleiotropy"))
  
  ### Return the dataset ###
  DATA$gene   = rownames(DATA)
  names(DATA) = c("Pleiotropy", "log_net_selection", "Pleiotropy_category", "gene")
  return(DATA)
}

### Build the eQTL phenotype dataset ###
build_eQTL_phenotype_dataset <- function( SNP_dataset, gene_dataset, eQTLs )
{
  ### Build the initial eQTL dataset ###
  DATA                  = data.frame(eQTLs)
  to_keep               = names(DATA)[!names(DATA)%in%c("Annotation","Putative_impact","gene")]
  DATA                  = DATA[,to_keep]
  DATA                  = merge(DATA, SNP_dataset, by="ID")
  DATA$net_selection    = gene_dataset[DATA$phenotype, "net_selection"]
  DATA$direct_selection = gene_dataset[DATA$phenotype, "direct_selection"]
  #DATA                  = drop_na(DATA)
  
  ### Calculate phenotype connectivity at the gene level ###
  Connectivity = c()
  for (p in unique(eQTLs$phenotype))
  {
    val          = length(unique(filter(eQTLs, phenotype==p)$gene))
    Connectivity = c(Connectivity, val)
  }
  names(Connectivity) = unique(eQTLs$phenotype)
  
  ### Calculate additional statistics ###
  DATA$Connectivity         = Connectivity[DATA$phenotype]
  DATA                      = DATA[!duplicated(DATA$phenotype),]
  DATA$log_net_selection    = log10(abs(DATA$net_selection))
  DATA$log_direct_selection = log10(abs(DATA$direct_selection))
  
  ### Build connectivity classes ###
  DATA$Connectivity_category = DATA$Connectivity
  DATA$Connectivity_category[DATA$Connectivity_category>1 & DATA$Connectivity_category<5] = 2
  DATA$Connectivity_category[DATA$Connectivity_category>=5] = 3
  DATA$Connectivity_category = as.character(DATA$Connectivity_category)
  DATA$Connectivity_category[DATA$Connectivity_category=="1"] = "Single eQTL"
  DATA$Connectivity_category[DATA$Connectivity_category=="2"] = "Low connectivity"
  DATA$Connectivity_category[DATA$Connectivity_category=="3"] = "High connectivity"
  DATA$Connectivity_category = factor(DATA$Connectivity_category, levels=c("Single eQTL","Low connectivity","High connectivity"))

  ### Build connectivity classes ###
  # DATA$Connectivity_category = DATA$Connectivity
  # DATA$Connectivity_category[DATA$Connectivity_category>1] = 2
  # DATA$Connectivity_category = as.character(DATA$Connectivity_category)
  # DATA$Connectivity_category[DATA$Connectivity_category=="1"] = "Single eQTL"
  # DATA$Connectivity_category[DATA$Connectivity_category=="2"] = "At least two eQTLs"
  # DATA$Connectivity_category = factor(DATA$Connectivity_category, levels=c("Single eQTL","At least two eQTLs"))
  # 
  ### Build parallel classes ###
  DATA$Parallel_category = as.character(DATA$Parallelism)
  DATA$Parallel_category[DATA$Parallel_category=="0"] = "No significant shift"
  DATA$Parallel_category[DATA$Parallel_category=="1"] = "Single shift"
  DATA$Parallel_category[DATA$Parallel_category=="2"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="3"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="4"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="5"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="6"] = "Partial or full parallelism"
  DATA$Parallel_category[DATA$Parallel_category=="7"] = "Partial or full parallelism"
  DATA$Parallel_category = factor(DATA$Parallel_category, levels=c("No significant shift","Single shift","Partial or full parallelism"))
  
  ### Return the dataset ###
  rownames(DATA) = DATA$phenotype
  return(DATA)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

SIZE_LIMIT = 10

#-----------------------------------#
# 1) Load and build input datasets  #
#-----------------------------------#
eQTLs                  = load_eQTLs("EXPRESSION")
SNP_dataset            = build_SNP_dataset()
gene_dataset           = build_gene_dataset(SNP_dataset, eQTLs)
eQTL_dataset           = build_eQTL_dataset(SNP_dataset, gene_dataset, eQTLs)
eQTL_carrier_dataset   = build_eQTL_carrier_dataset(eQTL_dataset, gene_dataset)
eQTL_phenotype_dataset = build_eQTL_phenotype_dataset(SNP_dataset, gene_dataset, eQTLs)
haplotype_blocks       = readRDS("./data/tribolium_haplotype_blocks/haplotype_blocks_CT_HD_G1_Tcas3.30.rds")
haplotype_blocks       = haplotype_blocks[,c("ID", "haplotype_block")]

#-----------------------------------#
# 2) Build the merged dataset       #
#-----------------------------------#
### Merge information ###
MERGED                         = merge(SNP_dataset, haplotype_blocks, by="ID")
MERGED$eQTL_carrier            = gene_dataset[MERGED$gene,"eQTL_carrier"]
MERGED$eQTL_carrier_category   = gene_dataset[MERGED$gene,"eQTL_carrier_category"]
MERGED$eQTL_phenotype          = gene_dataset[MERGED$gene,"eQTL_phenotype"]
MERGED$eQTL_phenotype_category = gene_dataset[MERGED$gene,"eQTL_phenotype_category"]
MERGED                         = merge(MERGED, eQTL_dataset[,c("ID", "beta", "phenotype", "Pleiotropy", "Pleiotropy_category")], by="ID", all.x=T)
### Update pleiotropy variables ###
MERGED$Pleiotropy[is.na(MERGED$Pleiotropy)]                   = 0
MERGED$Pleiotropy_category                                    = as.character(MERGED$Pleiotropy_category)
MERGED$Pleiotropy_category[is.na(MERGED$Pleiotropy_category)] = "No phenotype"
MERGED$Pleiotropy_category = factor(MERGED$Pleiotropy_category, levels=c("No phenotype", "No pleiotropy","Pleiotropy"))

#-----------------------------------#
# 3) Remove haplotype blocks with a #
#    size smaller than the limit    #
#-----------------------------------#
hb_size = table(MERGED$haplotype_block)
to_keep = names(hb_size[hb_size>=SIZE_LIMIT])
MERGED  = MERGED[MERGED$haplotype_block%in%to_keep,]

#-----------------------------------#
# 4) Save datasets in RDS format    #
#-----------------------------------#
saveRDS(SNP_dataset, file="./analyses/genetic_map_analysis/data/SNP_dataset.rds")
write.table(SNP_dataset, file="./analyses/genetic_map_analysis/data/SNP_dataset.csv", sep=";", col.names=T, row.names=F, quote=F)

saveRDS(gene_dataset, file="./analyses/genetic_map_analysis/data/gene_dataset.rds")
write.table(gene_dataset, file="./analyses/genetic_map_analysis/data/gene_dataset.csv", sep=";", col.names=T, row.names=F, quote=F)

saveRDS(MERGED, file="./analyses/genetic_map_analysis/data/merged_dataset.rds")
write.table(MERGED, file="./analyses/genetic_map_analysis/data/merged_dataset.csv", sep=";", col.names=T, row.names=F, quote=F)

