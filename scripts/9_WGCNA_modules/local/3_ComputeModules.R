#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_ComputeModules.R
# ------------------
# Find modules in a gene expression dataset with the Rpackage WGCNA.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")
library("WGCNA")

### Get the maximum absolute fitness correlation and p-value for a given power ###
get_fitness_correlation <- function( expression, fitness, power_value )
{
  #----------------------------------------------------------#
  # 1) Constructing the gene network and identifying modules #
  #----------------------------------------------------------#
  net = blockwiseModules(expression, power=power_value, TOMType="signed", minModuleSize=20,
                         reassignThreshold=0, mergeCutHeight=0.1, numericLabels=T,
                         pamRespectsDendro=F, saveTOMs=F, networkType = "signed",
                         saveTOMFileBase=paste0("./data/tribolium_modules/save"),
                         verbose=3, maxBlockSize=dim(expression)[2], deepSplit=4, randomSeed=1234)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs          = net$MEs
  geneTree     = net$dendrograms[[1]]
  #----------------------------------------------------------#
  # 2) Define numbers of genes and samples                   #
  #----------------------------------------------------------#
  nGenes   = ncol(expression)
  nSamples = nrow(expression)
  #----------------------------------------------------------#
  # 3) Recalculate MEs with color labels                     #
  #----------------------------------------------------------#
  MEs0              = moduleEigengenes(expression, moduleColors)$eigengenes
  MEs               = orderMEs(MEs0)
  moduleTraitCor    = cor(MEs, fitness, use="p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  index     = which(moduleTraitPvalue==min(moduleTraitPvalue))
  best_cor  = moduleTraitCor[index]
  best_pval = moduleTraitPvalue[index]
  best_nb   = length(unique(net$colors))
  return(c(power_value, best_cor, best_pval, best_nb))
}

### Get the best network and save it ###
get_and_save_best_network <- function( population, version, phenotype, expression, fitness, power_value )
{
  net = blockwiseModules(expression, power=power_value, TOMType="signed", minModuleSize=20,
                         reassignThreshold=0, mergeCutHeight=0.1, numericLabels=T,
                         pamRespectsDendro=F, saveTOMs=F, networkType = "signed",
                         saveTOMFileBase=paste0("./data/tribolium_modules/",population,"_",version,"_",phenotype),
                         verbose=3, maxBlockSize=dim(expression)[2], deepSplit=4, randomSeed=1234)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs          = net$MEs
  geneTree     = net$dendrograms[[1]]
  saveRDS(net, file=paste0("./data/tribolium_modules/network_",population,"_",version,"_",phenotype,".rds"))
  return(net)
}

### Extract significant modules ###
extract_significant_modules <- function( expression, fitness, net )
{
  #-----------------------------------------#
  # 1) Extract pertinent information        #
  #-----------------------------------------#
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs          = net$MEs
  geneTree     = net$dendrograms[[1]]
  #-----------------------------------------#
  # 2) Define numbers of genes and samples  #
  #-----------------------------------------#
  nGenes   = ncol(expression)
  nSamples = nrow(expression)
  #-----------------------------------------#
  # 3) Recalculate MEs with color labels    #
  #-----------------------------------------#
  MEs0              = moduleEigengenes(expression, moduleColors)$eigengenes
  MEs               = orderMEs(MEs0)
  moduleTraitCor    = cor(MEs, fitness, use="p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  #-----------------------------------------#
  # 4) Extract significant module names     #
  #-----------------------------------------#
  signif_index        = which(moduleTraitPvalue<0.05)
  signif_module_names = rownames(moduleTraitPvalue)[signif_index]
  signif_module_names = str_replace(signif_module_names, "ME", "")
  #-----------------------------------------#
  # 5) Extract the list of associated genes #
  #-----------------------------------------#
  color_list          = labels2colors(net$colors)
  names(color_list)   = names(net$colors)
  significant_modules = color_list[color_list%in%signif_module_names]
  significant_modules = data.frame(names(significant_modules), significant_modules)
  colnames(significant_modules) = c("gene","module")
  return(significant_modules)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation")

#-------------------------------------#
# 1) Read command line arguments      #
#-------------------------------------#
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)<4)
# {
#   stop("Please provide all command line arguments. Exit.", call.=FALSE)
# }
# POPULATION = args[1]
# VERSION    = args[2]
# PHENOTYPE  = args[3]
# POWER      = as.numeric(args[4])

POPULATION = "HD_G1"
VERSION    = "Tcas3.30"
PHENOTYPE  = "NOISE"
POWER      = 6

#-------------------------------------#
# 2) Load the dataset and the fitness #
#-------------------------------------#

### Load expression data ###
expr_filename        = paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_",PHENOTYPE,".txt")
expression           = read.table(expr_filename, h=T, sep="\t", check.names=F)
rownames(expression) = expression[,1]
expression           = t(expression[,-1])
gene_names           = colnames(expression)

### Load fitness data ###
fitness_filename  = paste0("./data/tribolium_phenotypes/Tribolium_castaneum_",POPULATION,"_",VERSION,"_FITNESS.txt")
fitness           = read.table(fitness_filename, h=T, sep="\t", check.names=F)
sample_names      = fitness[,1]
fitness           = data.frame(fitness$fitness)
rownames(fitness) = sample_names
colnames(fitness) = c("fitness")

#-------------------------------------#
# 3) Get and save the best network    #
#-------------------------------------#
BEST_NETWORK = get_and_save_best_network(POPULATION, VERSION, PHENOTYPE, expression, fitness, POWER)

#-------------------------------------#
# 4) Save significant modules         #
#-------------------------------------#
significant_modules = extract_significant_modules(expression, fitness, BEST_NETWORK)
write.table(significant_modules, file=paste0("./data/tribolium_modules/",POPULATION,"_",VERSION,"_",PHENOTYPE,"_significant_modules.txt"), col.names=T, row.names=F, sep="\t", quote=F)

