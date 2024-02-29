#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_CMH_tests.R
# -------------
# Run CMH tests on a set of SNPs.
# (HPC SCRIPT --> array wrapper)
#***************************************************************************

rm(list=ls())

library("tidyverse")

source("/scratch/project_2003847/lawstat-master/R/cmh.test.R")

### Run a CMH test for a given SNP ###
run_cochran_mantel_haenszel_test <- function( snp_id, data )
{
  X       = array(0, dim=c(2,2,7))
  x1      = select(filter(data, ID==snp_id & LINE=="L1"), c("AC", "NCALLED"))
  x2      = select(filter(data, ID==snp_id & LINE=="L2"), c("AC", "NCALLED"))
  x3      = select(filter(data, ID==snp_id & LINE=="L3"), c("AC", "NCALLED"))
  x4      = select(filter(data, ID==snp_id & LINE=="L5"), c("AC", "NCALLED"))
  x5      = select(filter(data, ID==snp_id & LINE=="L6"), c("AC", "NCALLED"))
  x6      = select(filter(data, ID==snp_id & LINE=="Mx1"), c("AC", "NCALLED"))
  x7      = select(filter(data, ID==snp_id & LINE=="Mx2"), c("AC", "NCALLED"))
  X[,,1]  = as.matrix(x1)
  X[,,2]  = as.matrix(x2)
  X[,,3]  = as.matrix(x3)
  X[,,4]  = as.matrix(x4)
  X[,,5]  = as.matrix(x5)
  X[,,6]  = as.matrix(x6)
  X[,,7]  = as.matrix(x7)
  test    = cmh.test(X)
  return(as.vector(c(test$parameter["Pooled Odd Ratio"][[1]], test$parameter["p-value"][[1]])))
}


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
SIZE  = as.numeric(args[1])
INDEX = as.numeric(args[2])
PATH  = args[3]

print(paste0("SIZE : ",SIZE))
print(paste0("INDEX : ",INDEX))
print(paste0("PATH : ",PATH))

setwd(PATH)

#--------------------------------#
# 2) Load and prepare datasets   #
#--------------------------------#
DATA          = readRDS("/scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_afc/AF.rds")
DISTINCT_DATA = distinct(DATA, ID, .keep_all=T)
DISTINCT_DATA = select(DISTINCT_DATA, c("ID", "POS", "CHROM"))
START         = (INDEX-1)*SIZE+1
END           = INDEX*SIZE
if(END > dim(DISTINCT_DATA)[1])
{
  END = dim(DISTINCT_DATA)[1]
}

#--------------------------------#
# 3) Run the tests               #
#--------------------------------#
res        = lapply(DISTINCT_DATA$ID[START:END], function(x) {run_cochran_mantel_haenszel_test(x, DATA)})
res        = as.data.frame(do.call(rbind, res))
names(res) = c("OR", "PVALUE")
res$ID     = DISTINCT_DATA$ID[START:END]
res$POS    = DISTINCT_DATA$POS[START:END]
res$CHROM  = DISTINCT_DATA$CHROM[START:END]

#--------------------------------#
# 4) Save the result             #
#--------------------------------#
filename = paste0(INDEX,"_cmh_tests.txt")
write.table(res, file=filename, col.names=T, row.names=F, quote=F)
system(paste0("cp ",filename," /scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_afc/cmh_tests/",filename))

