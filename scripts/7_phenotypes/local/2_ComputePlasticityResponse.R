#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_ComputePlasticityResponse.R
# -----------------------------
# Compute plasticity response between CT and HD.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("edgeR")
library("limma")


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#----------------------------------------------#
# 1) Read command line arguments               #
#----------------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
POPULATION = args[1]
VERSION    = args[2]

# POPULATION = "CT_HD_G1"
# VERSION    = "Tcas3.30"

#----------------------------------------------#
# 2) Loading samples                           #
#----------------------------------------------#
samples                = read.table(paste0("./data/tribolium_bam/samples_",POPULATION,"_",VERSION,".csv"), h=T, sep=";")
samples$fem            = as.factor(samples$fem)
samples$line           = as.factor(samples$line)
samples$fullsib_family = as.factor(samples$fullsib_family)
samples$halfsib_family = as.factor(samples$halfsib_family)
samples$source_env     = as.factor(samples$source_env)
samples$target_env     = as.factor(samples$target_env)
samples$run_index      = as.factor(samples$run_index)
samples$batch_index    = as.factor(samples$batch_index)

#----------------------------------------------#
# 3) Loading read counts                       #
#----------------------------------------------#
d       = read.table(paste0("./data/tribolium_counts/Tribolium_castaneum_",POPULATION,"_",VERSION,"_read_counts_READY.txt"), sep="\t", h=T, check.names=F)
gene_id = d[,1]
X       = d[,samples$sample]

#----------------------------------------------#
# 4) Calculate expression differential between #
#    each HD individual and its mean CT family #
#----------------------------------------------#
HD_samples = samples[samples$target_env=="HD",]
dX         = X[,HD_samples$sample]
for(i in seq(1,dim(HD_samples)[1]))
{
  sample      = HD_samples$sample[i]
  family      = HD_samples$fem[i]
  CT_samples  = samples[samples$target_env=="CT" & samples$fem==family,]
  mean_X      = rowMeans(X[,CT_samples$sample])
  dX[,sample] = dX[,sample]-mean_X
}

#----------------------------------------------#
# 5) Quantile normalization                    #
#----------------------------------------------#
for(i in 1:dim(dX)[1])
{
  mat    = dX[i,]
  mat    = rank(mat, ties.method="average")
  mat    = qnorm(mat/(ncol(dX)+1))
  dX[i,] = mat
}
rm(i, mat)
dX = as.data.frame(dX)
dX = cbind(gene_id, dX)

#----------------------------------------------#
# 6) Save the normalized and adjusted dataset  #
#----------------------------------------------#
write.table(dX, file=paste0("./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_PLASTICITY.txt"), quote=F, sep="\t", row.names=F, col.names=T)
system(paste0("gzip -c ./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_PLASTICITY.txt > ./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_PLASTICITY.txt.gz"))

