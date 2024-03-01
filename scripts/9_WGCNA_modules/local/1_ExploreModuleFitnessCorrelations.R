#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_ExploreModuleFitnessCorrelations.R
# ------------------------------------
# Use Rpackage WGCNA to explore functional modules of gene expression and
# their correlation to fitness.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")
library("WGCNA")

### Get the maximum absolute fitness correlation and p-value for a given power ###
get_fitness_correlation <- function( expression, fitness, power_value, cor_method )
{
  #----------------------------------------------------------#
  # 1) Constructing the gene network and identifying modules #
  #----------------------------------------------------------#
  net = blockwiseModules(expression, power=power_value, TOMType="signed", minModuleSize=20,
                         reassignThreshold=0, mergeCutHeight=0.1, numericLabels=T,
                         corType=cor_method,
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
get_and_save_best_network <- function( population, version, phenotype, expression, fitness, power_value, cor_method )
{
  net = blockwiseModules(expression, power=power_value, TOMType="signed", minModuleSize=20,
                         corType=cor_method,
                         reassignThreshold=0, mergeCutHeight=0.1, numericLabels=T,
                         pamRespectsDendro=F, saveTOMs=T, networkType = "signed",
                         saveTOMFileBase=paste0("./data/tribolium_modules/",population,"_",version,"_",phenotype),
                         verbose=3, maxBlockSize=dim(expression)[2], deepSplit=4, randomSeed=1234)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs          = net$MEs
  geneTree     = net$dendrograms[[1]]
  save(MEs, moduleLabels, moduleColors, geneTree, file=paste0("./data/tribolium_modules/modules_",population,"_",version,"_",phenotype,"_fitness_raw.RData"))
  return(net)
}

### Plot the correlation matrix with fitness ###
plot_fitness_correlation <- function( expression, fitness, net )
{
  #----------------------------------------------------------#
  # 1) Extract pertinent information                         #
  #----------------------------------------------------------#
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
  #----------------------------------------------------------#
  # 4) Make the plot                                         #
  #----------------------------------------------------------#
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(fitness),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}

### Plot the repartition of genes per module ###
plot_nb_genes_per_module <- function( net )
{
  TB = table(labels2colors(net$colors))
  TB = sort(TB, decreasing=T)
  barplot(TB, las=2, col=names(TB), main="Repartition of genes in modules")
}


##################
#      MAIN      #
##################

#-------------------------------------#
# 1) Read command line arguments      #
#-------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
COR_METHOD      = args[4]
PHENOTYPE       = args[5]
setwd(REPOSITORY_PATH)

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
# 3) Explore fitness correlations     #
#-------------------------------------#
fitness_cor_data = c()
for (power_value in 1:30)
{
  res = get_fitness_correlation(expression, fitness, power_value, COR_METHOD)
  print(paste(power_value, res))
  fitness_cor_data = rbind(fitness_cor_data, res)
}
fitness_cor_data        = as.data.frame(fitness_cor_data)
names(fitness_cor_data) = c("power", "best_cor", "best_pval", "nb_modules")

#-------------------------------------#
# 4) Compute the scale-free fit       #
#-------------------------------------#
powers               = c(1:30)
sft                  = pickSoftThreshold(expression, powerVector=powers, verbose=5)
fitness_cor_data$fit = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]

#-------------------------------------#
# 5) Save the exploration result      #
#-------------------------------------#
write.table(fitness_cor_data, file=paste0("./data/tribolium_modules/power_exploration_",POPULATION,"_",VERSION,"_",COR_METHOD,"_",PHENOTYPE,".txt"), quote=F, sep="\t", row.names=F, col.names=T)

#-------------------------------------#
# 6) Plot the result                  #
#-------------------------------------#
fitness_cor_data = read.table(paste0("./data/tribolium_modules/power_exploration_",POPULATION,"_",VERSION,"_",COR_METHOD,"_",PHENOTYPE,".txt"), sep="\t", h=T)
THRESHOLD = 6

par(mfrow=c(2,2))
plot(fitness_cor_data$power, fitness_cor_data$best_cor, pch=20, xlab="Soft threshold (power)", ylab="Correlation", main="Module with best correlation to fitness")
lines(fitness_cor_data$power, fitness_cor_data$best_cor)
abline(v=THRESHOLD, col="red", lty=2)

plot(fitness_cor_data$power, fitness_cor_data$best_pval, pch=20, xlab="Soft threshold (power)", ylab="P-value", main="Associated p-value", log="y")
lines(fitness_cor_data$power, fitness_cor_data$best_pval)
abline(v=THRESHOLD, col="red", lty=2)
abline(h=0.05, col="purple")
plot(fitness_cor_data$power, fitness_cor_data$nb_modules, pch=20, xlab="Soft threshold (power)", ylab="Number of modules", main="Number of modules found")
lines(fitness_cor_data$power, fitness_cor_data$nb_modules)
abline(v=THRESHOLD, col="red", lty=2)

plot(fitness_cor_data$power, fitness_cor_data$fit, pch=20, xlab="Soft threshold (power)", ylab="Scale-free model fit", main="Scale-free model fitting")
lines(fitness_cor_data$power, fitness_cor_data$fit)
abline(v=THRESHOLD, col="red", lty=2)

