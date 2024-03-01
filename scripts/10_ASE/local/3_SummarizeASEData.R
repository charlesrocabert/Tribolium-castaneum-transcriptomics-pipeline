#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
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

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
SUFFIX          = args[2]
VERSION         = args[3]

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

