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
  AFC_HD = AFC_HD[,c("ID", "POS")]#, "L1_AFC", "L2_AFC", "L3_AFC", "L5_AFC", "L6_AFC", "Mx1_AFC", "Mx2_AFC", "L1_AFC_abs", "L2_AFC_abs", "L3_AFC_abs", "L5_AFC_abs", "L6_AFC_abs", "Mx1_AFC_abs", "Mx2_AFC_abs")]
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
  ANNOTATION_SELECT = select(ANNOTATION, c("ID", "CHROM", "Annotation", "Putative_impact", "gene"))
  
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
  DATA$L1_shift  = as.numeric(DATA$ID%in%DATA$ID[L1_SIGNIF$index])
  DATA$L2_shift  = as.numeric(DATA$ID%in%DATA$ID[L2_SIGNIF$index])
  DATA$L3_shift  = as.numeric(DATA$ID%in%DATA$ID[L3_SIGNIF$index])
  DATA$L5_shift  = as.numeric(DATA$ID%in%DATA$ID[L5_SIGNIF$index])
  DATA$L6_shift  = as.numeric(DATA$ID%in%DATA$ID[L6_SIGNIF$index])
  DATA$Mx1_shift = as.numeric(DATA$ID%in%DATA$ID[Mx1_SIGNIF$index])
  DATA$Mx2_shift = as.numeric(DATA$ID%in%DATA$ID[Mx2_SIGNIF$index])
  
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
  #DATA = merge(DATA, PHENOTYPE_POSITIONS, by="gene")
  
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


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

SIZE_LIMIT = 1

#-----------------------------------#
# 1) Load and build input datasets  #
#-----------------------------------#
eQTLs            = load_eQTLs("EXPRESSION")
SNP_dataset      = build_SNP_dataset()
haplotype_blocks = readRDS("./data/tribolium_haplotype_blocks/haplotype_blocks_CT_HD_G1_Tcas3.30.rds")
haplotype_blocks = haplotype_blocks[,c("ID", "haplotype_block")]

#-----------------------------------#
# 2) Build the merged dataset       #
#-----------------------------------#
HB_dataset                = merge(SNP_dataset, haplotype_blocks, by="ID")
hb_size                   = table(HB_dataset$haplotype_block)
HB_dataset$hb_size        = hb_size[HB_dataset$haplotype_block]
HB_dataset$eQTL_carrier   = as.numeric(HB_dataset$gene %in% eQTLs$gene)
HB_dataset$eQTL_phenotype = as.numeric(HB_dataset$gene %in% eQTLs$phenotype)

#-----------------------------------#
# 3) Remove haplotype blocks with a #
#    size smaller than the limit    #
#-----------------------------------#
HB_dataset = filter(HB_dataset, hb_size >= SIZE_LIMIT)

#-----------------------------------#
# 4) Save datasets in RDS format    #
#-----------------------------------#
saveRDS(HB_dataset, file="./analyses/genetic_map_analysis/data/haplotype_block_dataset.rds")
write.table(HB_dataset, file="./analyses/genetic_map_analysis/data/haplotype_block_dataset.csv", sep=";", col.names=T, row.names=F, quote=F)

length(unique(HB_dataset$haplotype_block))
4921/50
99*50
