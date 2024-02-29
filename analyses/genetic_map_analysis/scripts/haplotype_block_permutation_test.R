#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# haplotype_block_permutation_test.R
# ----------------------------------
# Run a permutation test on a given haplotype block.
# (HPC SCRIPT --> array wrapper)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Aggregator function ###
aggregator_function <- function( x, aggregator )
{
  if (aggregator == "sum")
  {
    return(sum(x, na.rm=T))
  }
  if (aggregator == "mean")
  {
    return(mean(x, na.rm=T))
  }
  if (aggregator == "median")
  {
    return(median(x, na.rm=T))
  }
}

### Run one permutations test ###
run_permutations_test <- function( dataset, target_variable, source_variable, source_list, nb_reps, aggregator )
{
  x              = dataset[,target_variable]
  pos            = which(dataset[,source_variable]%in%source_list)
  original_score = aggregator_function(x[pos], aggregator)
  permutations   = t(sapply(1:nb_reps,function(i){
    score = aggregator_function(sample(x)[pos], aggregator)
    return(score)
  }))
  qval = sum(permutations<original_score)/nb_reps
  return(1-qval)
}


##################
#      MAIN      #
##################

#--------------------------------#
# 1) Read command line arguments #
#--------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
VARIABLE   = args[1]
REPS       = as.numeric(args[3])
AGGREGATOR = args[2]
SIZE       = as.numeric(args[4])
INDEX      = as.numeric(args[5])
PATH       = args[6]

print(paste0("VARIABLE   : ",VARIABLE))
print(paste0("AGGREGATOR : ",AGGREGATOR))
print(paste0("REPS       : ",REPS))
print(paste0("SIZE       : ",SIZE))
print(paste0("INDEX      : ",INDEX))
print(paste0("PATH       : ",PATH))

setwd(PATH)

#--------------------------------#
# 2) Load and prepare datasets   #
#--------------------------------#
DATA     = readRDS("/scratch/project_2003847/Tribolium_castaneum_haplotype_blocks/haplotype_block_dataset.rds")
HB_NAMES = unique(DATA$haplotype_block)
START    = (INDEX-1)*SIZE+1
END      = INDEX*SIZE
if(END > length(HB_NAMES))
{
  END = length(HB_NAMES)
}

#--------------------------------#
# 3) Run the permutation tests   #
#--------------------------------#
pvalues = c()
for (i in seq(START, END))
{
  hb_name = HB_NAMES[i]
  ID_list = filter(DATA, haplotype_block==hb_name)$ID
  pval    = run_permutations_test(DATA, VARIABLE, "ID", ID_list, REPS, AGGREGATOR)
  pvalues = c(pvalues, pval)
}
RES        = data.frame(seq(START, END), HB_NAMES[START:END], rep(VARIABLE, length(pvalues)), rep(REPS, length(pvalues)), rep(AGGREGATOR, length(pvalues)), pvalues)
names(RES) = c("index", "haplotype_block", "variable", "nb_rep", "aggregator", "pvalue")

#--------------------------------#
# 4) Save the result             #
#--------------------------------#
filename = paste0(INDEX, "_", VARIABLE, ".txt")
write.table(RES, file=filename, sep=";", col.names=T, row.names=F, quote=F)
system(paste0("cp ", filename, " /scratch/project_2003847/Tribolium_castaneum_haplotype_blocks/output/", filename))

