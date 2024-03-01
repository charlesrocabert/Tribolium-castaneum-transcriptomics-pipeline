#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 4_PlotSelectedModules.R
# -----------------------
# Plot the result of the modules detection.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")
library("WGCNA")

### Plot the correlation matrix with fitness ###
plot_fitness_correlation <- function( expression, fitness, net )
{
  #----------------------------------------#
  # 1) Extract pertinent information       #
  #----------------------------------------#
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs          = net$MEs
  geneTree     = net$dendrograms[[1]]
  gene_count   = table(labels2colors(net$colors))
  #----------------------------------------#
  # 2) Define numbers of genes and samples #
  #----------------------------------------#
  nGenes   = ncol(expression)
  nSamples = nrow(expression)
  #----------------------------------------#
  # 3) Recalculate MEs with color labels   #
  #----------------------------------------#
  MEs0              = moduleEigengenes(expression, moduleColors)$eigengenes
  MEs               = orderMEs(MEs0)
  raw_colors        = str_replace(names(MEs), "ME", "")
  new_names         = paste0(gene_count[raw_colors], " genes in ", raw_colors)
  moduleTraitCor    = cor(MEs, fitness, use="p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  index             = which(moduleTraitPvalue==min(moduleTraitPvalue))
  best_cor          = moduleTraitCor[index]
  best_pval         = moduleTraitPvalue[index]
  best_nb           = length(unique(net$colors))
  #----------------------------------------#
  # 4) Make the plot                       #
  #----------------------------------------#
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar=c(4,12,4,4))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(fitness),
                 yLabels = names(MEs),
                 ySymbols = new_names,
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-fitness relationships")
}
plot_fitness_correlation(expression, fitness, net)

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

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
PHENOTYPE       = args[4]

#--------------------------------#
# 2) Load the datasets           #
#--------------------------------#

### Load the network ###
net = readRDS(paste0("./data/tribolium_modules/network_",POPULATION,"_",VERSION,"_",PHENOTYPE,".rds"))

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

#--------------------------------#
# 3) Make plots                  #
#--------------------------------#
plot_fitness_correlation(expression, fitness, net)
#plot_nb_genes_per_module(net)

