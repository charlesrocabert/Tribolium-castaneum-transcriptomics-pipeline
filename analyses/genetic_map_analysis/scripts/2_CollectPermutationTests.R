#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_CollectPermutationTests.R
# ---------------------------
# Collect the result of permutation tests.
# (HPC SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")


##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
VARIABLE = args[1]
N        = as.numeric(args[2])

print(paste0("VARIABLE : ",VARIABLE))
print(paste0("N        : ",N))

#--------------------------------#
# 2) Import permutation results  #
#--------------------------------#
MERGED = data.frame()
for (i in seq(1,N))
{
    print(paste0("> Importing permutation test ", i, " of ", N))
    d = read.table(paste0("/scratch/project_2003847/Tribolium_castaneum_haplotype_blocks/output/", i, "_", VARIABLE, ".txt"), sep=";", header=T)
    MERGED = rbind(MERGED, d)
}

#--------------------------------#
# 3) Save the result             #
#--------------------------------#
write.table(MERGED, file=paste0("/scratch/project_2003847/Tribolium_castaneum_haplotype_blocks/", VARIABLE, "_merged.txt"), sep=";", row.names=F, col.names=T, quote=F)
