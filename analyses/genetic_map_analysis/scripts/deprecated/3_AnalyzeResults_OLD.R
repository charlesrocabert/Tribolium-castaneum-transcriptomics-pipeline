#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_MergeResults.R
# ----------------
# Merge locally all permutation results in a single file.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("mltools")
library("qqman")
library("ggpubr")
library("ggvenn")
library("RColorBrewer")
library("plotROC")

### Calculate the FDR ###
calculate_fdr <- function( p )
{
  fdr = p.adjust(p, method="fdr")
  return(fdr)
}

### Load permutation results and merge with additional data ###
load_permutation_results <- function( hb_data, classification_score )
{
  #---------------------------------------#
  # 1) Load and merge permutation results #
  #---------------------------------------#
  PREFIXES = c("isShifting", "isParallel", "evol_response_pool", "eQTL_carrier", "eQTL_phenotype")
  MERGED   = data.frame()
  first    = T
  for (prefix in PREFIXES)
  {
    d = read.table(paste0("./analyses/genetic_map_analysis/data/",prefix,"_merged.txt"), sep=";", h=T)
    if (first)
    {
      MERGED        = data.frame(d[,c("index", "haplotype_block", "pvalue")])
      MERGED$fdr    = calculate_fdr(MERGED$pvalue)
      if (classification_score == "pvalue")
      {
        MERGED$signif = as.numeric(MERGED$pvalue < 0.05)
      } else if (classification_score == "FDR") {
        MERGED$signif = as.numeric(MERGED$fdr < 0.05)
      }
      names(MERGED) = c("index", "haplotype_block", paste0(prefix,"_pvalue"), paste0(prefix,"_fdr"), prefix)
      first         = F
    } else {
      d     = d[,c("index", "pvalue")]
      d$fdr = calculate_fdr(d$pvalue)
      if (classification_score == "pvalue")
      {
        d$signif = as.numeric(d$pvalue < 0.05)
      } else if (classification_score == "FDR") {
        d$signif = as.numeric(d$fdr < 0.05)
      }
      names(d) = c("index", paste0(prefix,"_pvalue"), paste0(prefix,"_fdr"), prefix)
      MERGED   = merge(MERGED, d, by="index")
    }
  }
  #---------------------------------------#
  # 2) Collect additional information     #
  #---------------------------------------#
  ### 2.1) Chromosome location ###
  hb_chr            = hb_data[!duplicated(hb_data$haplotype_block),c("haplotype_block", "gene_chr")]
  row.names(hb_chr) = hb_chr$haplotype_block
  MERGED$chr        = hb_chr[as.character(MERGED$haplotype_block),"gene_chr"]
  MERGED            = MERGED[order(MERGED$chr, MERGED$haplotype_block),]
  ### 2.2) Custom index per chromosome ###
  current_chr = MERGED$chr[1]
  counter     = 1
  chr_index   = c()
  for(i in seq(1,dim(MERGED)[1]))
  {
    new_chr = MERGED$chr[i]
    if (current_chr == new_chr)
    {
      chr_index = c(chr_index, counter)
      counter   = counter+1
    } else {
      current_chr = new_chr
      counter     = 1
      chr_index = c(chr_index, counter)
    }
  }
  MERGED$chr_index = chr_index
  ### 2.2) Mean log(|Sg|) ###
  mean_sg = tapply(hb_data$log_net_selection, hb_data$haplotype_block, mean)
  MERGED$mean_sg = mean_sg[as.character(MERGED$haplotype_block)]
  ### 2.3) Max log(|Sg|) ###
  max_sg = tapply(hb_data$log_net_selection, hb_data$haplotype_block, max)
  MERGED$max_sg = max_sg[as.character(MERGED$haplotype_block)]
  ### 2.4) Max log(|Sg|) gene
  ### 2.5) Dominant gene
  ### 2.6) Mean abs(AFC) per line
  LINES = c("L1", "L2", "L3", "L5", "L6", "Mx1", "Mx2")
  for(line in LINES)
  {
    variable     = paste0(line,"_AFC_abs")
    new_variable = paste0("mean_",line,"_AFC_abs")
    mean_afc_abs = tapply(hb_data[,variable], hb_data$haplotype_block, mean)
    MERGED[,new_variable] = mean_afc_abs[as.character(MERGED$haplotype_block)]
  }
  ### 2.7) Max abs(AFC) per line
  LINES = c("L1", "L2", "L3", "L5", "L6", "Mx1", "Mx2")
  for(line in LINES)
  {
    variable     = paste0(line,"_AFC_abs")
    new_variable = paste0("max_",line,"_AFC_abs")
    mean_afc_abs = tapply(hb_data[,variable], hb_data$haplotype_block, max)
    MERGED[,new_variable] = mean_afc_abs[as.character(MERGED$haplotype_block)]
  }
  
  ### TEMPORARY ###
  #is_shifting = tapply(hb_data$isShifting, hb_data$haplotype_block, max)
  #MERGED$isShifting = is_shifting[as.character(MERGED$haplotype_block)]
  #is_parallel = tapply(hb_data$isParallel, hb_data$haplotype_block, max)
  #MERGED$isParallel = is_parallel[as.character(MERGED$haplotype_block)]
  
  shifting_ratio = tapply(hb_data$isShifting, hb_data$haplotype_block, mean)
  MERGED$shifting_ratio = shifting_ratio[as.character(MERGED$haplotype_block)]
  parallel_ratio = tapply(hb_data$isParallel, hb_data$haplotype_block, mean)
  MERGED$parallel_ratio = parallel_ratio[as.character(MERGED$haplotype_block)]
  DE_ratio = tapply(hb_data$evol_response_pool, hb_data$haplotype_block, mean)
  MERGED$DE_ratio = DE_ratio[as.character(MERGED$haplotype_block)]
  # is_eQTL = tapply(hb_data$eQTL_carrier, hb_data$haplotype_block, max)
  # MERGED$eQTL_carrier = is_eQTL[as.character(MERGED$haplotype_block)]
  # is_pheno = tapply(hb_data$eQTL_phenotype, hb_data$haplotype_block, max)
  # MERGED$eQTL_phenotype = is_pheno[as.character(MERGED$haplotype_block)]
  #---------------------------------------#
  # 3) Return the final merged data       #
  #---------------------------------------#
  return(MERGED)
}

### Run permutation test ###
run_permutation_test <- function( dataset, variable, index_list, nb_reps )
{
  x              = dataset[,variable]
  pos            = which(dataset$index%in%index_list)
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

HB_DATA   = readRDS(file="./analyses/genetic_map_analysis/data/merged_dataset.rds")
PERM_DATA = load_permutation_results(HB_DATA, "FDR")

# unique(HB_DATA$gene_chr)
# 
# hb_chr = HB_DATA[!duplicated(HB_DATA$haplotype_block),c("haplotype_block", "gene_chr")]
# row.names(hb_chr) = hb_chr$haplotype_block
# PERM_DATA$chr = hb_chr[as.character(PERM_DATA$haplotype_block),"gene_chr"]
# PERM_DATA$chr_int = as.numeric(as.factor(PERM_DATA$chr))
# PERM_DATA$isShifting_fdr = calculate_fdr(PERM_DATA$isShifting_pvalue)
# manhattan(PERM_DATA, chr="chr_int", bp="chr_index", snp="haplotype_block", p="isShifting_pvalue", ylim=c(0,10), chrlabs=unique(PERM_DATA$chr), las=2, xlab="")
# qq(PERM_DATA$isShifting_fdr)
# # fisher.test(table(PERM_DATA$isParallel, PERM_DATA$eQTL_phenotype))
# # 
# plot(-log10(calculate_fdr(PERM_DATA$isShifting_pvalue)), type="h")
# abline(h=-log10(0.05), col="red")
# PERM_DATA$isParallel_pvalue

#PERM_DATA$haplotype_block

# comp = list(c("0", "1"))
# ggplot(PERM_DATA, aes(as.factor(evol_response_pool), mean_Mx1_AFC_abs)) +
#   geom_boxplot() +
#   stat_compare_means(comparisons=comp, method="wilcox.test", label="p.signif") +
#   theme_classic()

# SHIFTING = filter(PERM_DATA, isShifting==1)$index
# DE = filter(PERM_DATA, evol_response_pool==1)$index
# EQTL = filter(PERM_DATA, eQTL_carrier==1)$index
# PHENO = filter(PERM_DATA, eQTL_phenotype==1)$index
# 
# a = list("Shifting"  = as.character(filter(PERM_DATA, isShifting==1)$index),
#          #"Parallel"  = as.character(filter(PERM_DATA, isParallel==1)$index),
#          "DE"        = as.character(filter(PERM_DATA, evol_response_pool==1)$index),
#          "eQTL"      = as.character(filter(PERM_DATA, eQTL_carrier==1)$index),
#          "Regulated" = as.character(filter(PERM_DATA, eQTL_phenotype==1)$index))
# ggvenn(a)
# a
# index_list = filter(PERM_DATA, isParallel==1)$index


PREFIXES = c("isShifting", "isParallel", "evol_response_pool", "eQTL_carrier", "eQTL_phenotype")
N        = length(PREFIXES)
for(i in seq(1,N-1))
{
  for(j in seq(i+1,N))
  {
    pval1 = fisher.test(table(PERM_DATA[,PREFIXES[i]], PERM_DATA[,PREFIXES[j]]))$p.value
    #pval2 = run_permutation_test(PERM_DATA, PREFIXES[i], PERM_DATA[PERM_DATA[,PREFIXES[j]]==1,"index"], 1000)
    print(paste(PREFIXES[i], PREFIXES[j], pval1, "(Fisher T)"))
    #print(paste(PREFIXES[i], PREFIXES[j], pval2, "(Perm)"))
    print(cor.test(PERM_DATA[,PREFIXES[i]], PERM_DATA[,PREFIXES[j]])$p.value)
    #print("-----------------------")
  }
}
# plot(PERM_DATA$shifting_ratio, PERM_DATA$DE_ratio, log="xy")
# plot(PERM_DATA$DE_ratio)
# cor.test(log10(PERM_DATA$shifting_ratio), log10(PERM_DATA$DE_ratio), method="spearman")
# table(PERM_DATA$isParallel)
# table(PERM_DATA$evol_response_pool)


# summary(lm(evol_response_pool~isParallel+eQTL_carrier, data=PERM_DATA))
# 
# summary(lm(DE_ratio~parallel_ratio*shifting_ratio, data=PERM_DATA))
# 
# p = ggplot(PERM_DATA, aes(d=isParallel, m=DE_ratio)) + geom_roc()
# p
# calc_auc(p)$AUC


# model_glm = glm(evol_response_pool~isParallel+eQTL_carrier, family="binomial", data=PERM_DATA)
# summary(model_glm)

# Predictions on the training set

# predictTrain = predict(model_glm, data = PERM_DATA, type = "response")


# Confusion matrix on training data

# table(PERM_DATA$evol_response_pool, predictTrain >= 0.5)
# cor.test(PERM_DATA$evol_response_pool, PERM_DATA$isParallel, method="spearman")
# table(PERM_DATA$evol_response_pool, PERM_DATA$isParallel)
# chisq.test(PERM_DATA$evol_response_pool, PERM_DATA$isShifting)
# 
# (114+268)/nrow(train) #Accuracy - 91%


#Predictions on the test set

# predictTest = predict(model_glm, newdata = test, type = "response")
# 

# Confusion matrix on test set

# table(test$approval_status, predictTest >= 0.5)
# 
# 158/nrow(test) #Accuracy - 88%

filter(PERM_DATA, isParallel==1)
