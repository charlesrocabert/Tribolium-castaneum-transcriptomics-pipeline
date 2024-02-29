#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_ComputePhenotypicNoise.R
# --------------------------
# Estimate phenotypic noise (Ve) in HD.
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

# POPULATION = "HD_G1"
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
d           = read.table(paste0("./data/tribolium_counts/Tribolium_castaneum_",POPULATION,"_",VERSION,"_read_counts_READY.txt"), sep="\t", h=T, check.names=F)
gene_id     = d[,1]
X           = d[,samples$sample]
rownames(X) = gene_id

#---------------------------------------------#
# 4) For each gene:                           #
#    - Remove family structure effects,       #
#    - Compute variances per family.          #
#---------------------------------------------#
X_residuals = c()
X_noise     = c()
for(i in 1:dim(X)[1])
{
  if (i%%1000==0)
  {
    writeLines(paste0(">>> Dealing with gene ",gene_id[i], " (", round(i/dim(X)[1]*100, 2), "%)"))
  }
  model       = lm(t(X[i,])~fullsib_family+halfsib_family, data=samples)
  D           = data.frame(t(X[i,]), predict(model), residuals(model))
  names(D)    = c("observed", "predicted", "residuals")
  D           = D[samples$sample,]
  D$sample    = samples$sample
  D$family    = samples$fem
  fam_sd      = tapply(D$residuals, D$family, sd)
  fam_mean    = tapply(D$residuals, D$family, mean)
  D$cv        = fam_sd[D$family]/fam_mean[D$family]
  X_residuals = rbind(X_residuals, D$residuals)
  X_noise     = rbind(X_noise, D$cv)
}
colnames(X_residuals) = samples$sample
rownames(X_residuals) = gene_id
colnames(X_noise)     = samples$sample
rownames(X_noise)     = gene_id

#---------------------------------------------#
# 5) Quantile normalization                   #
#---------------------------------------------#
for(i in 1:dim(X_noise)[1])
{
  mat         = X_noise[i,]
  mat         = rank(mat, ties.method="average")
  mat         = qnorm(mat/(ncol(X_noise)+1))
  X_noise[i,] = mat
}
rm(i, mat)
X_noise = as.data.frame(X_noise)
X_noise = cbind(gene_id, X_noise)

#---------------------------------------------#
# 6) Save the normalized and adjusted dataset #
#---------------------------------------------#
write.table(X_noise, file=paste0("./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_NOISE.txt"), quote=F, sep="\t", row.names=F, col.names=T)
system(paste0("gzip -c ./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_NOISE.txt > ./data/tribolium_phenotypes/Tribolium_castaneum_HD_G1_",VERSION,"_NOISE.txt.gz"))

