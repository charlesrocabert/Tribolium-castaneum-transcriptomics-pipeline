#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_SummarizeASEData.R
# --------------------
# Summarize ASE data at the population level and save the result.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
SUFFIX  = args[1]
VERSION = args[2]

# SUFFIX  = "HD_G1"
# VERSION = "Tcas3.30"

#--------------------------------#
# 2) Load the ASE datasets       #
#--------------------------------#
het_snps = readRDS(paste0("./data/tribolium_ASE/",SUFFIX,"_",VERSION,"_HET_SNPs.rds"))
sig_snps = readRDS(paste0("./data/tribolium_ASE/",SUFFIX,"_",VERSION,"_SIGNIF_SNPs.rds"))

#--------------------------------#
# 3) Summarize datasets          #
#--------------------------------#
X     = rowSums(het_snps[,2:M1])
Y     = rowSums(sig_snps[,2:M2])
Xn    = rowSums(het_snps[,2:M1])/M1
Yn    = rowSums(sig_snps[,2:M2])/M2
Xecdf = ecdf(X)(X)
Yecdf = ecdf(Y)(Y)

SUMMARY = data.frame(het_snps$snp_id, X, Xn, Xecdf, Y, Yn, Yecdf)
names(SUMMARY) = c("ID", "Het_count", "Het_count_norm", "Het_count_ecdf", "ASE_count", "ASE_count_norm", "ASE_count_ecdf")

write.table(SUMMARY, file=paste0("./data/tribolium_ASE/",SUFFIX,"_",VERSION,"_ASE_SNP_table.txt"), col.names=T, row.names=F, quote=F, sep="\t")

# par(mfrow=c(1,2))
# M1 = dim(het_snps)[2]
# hist(rowSums(het_snps[,2:M1]), main="Distribution of het count\namong SNPs", xlab="Nb heterozygote individuals")
# abline(v=130, col="red")
# M2 = dim(sig_snps)[2]
# hist(rowSums(sig_snps[,2:M2]), main="Distribution of signif count\namong SNPs", xlab="Nb significant individuals", freq=F)
