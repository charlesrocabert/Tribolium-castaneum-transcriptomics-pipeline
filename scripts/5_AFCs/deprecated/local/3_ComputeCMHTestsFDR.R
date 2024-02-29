#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_ComputeCMHTestsFDR.R
# ----------------------
# Compute FDR on all CMH tests.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

### Calculate the FDR ###
calculate_fdr <- function( p )
{
  fdr = p.adjust(p, method="fdr")
  return(fdr)
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#--------------------------#
# 1) Load the test dataset #
#--------------------------#
CMH_TESTS        = read.table("./data/tribolium_afc/CMH_TESTS.txt", h=T, sep=" ")
CMH_TESTS$PVALUE = as.numeric(CMH_TESTS$PVALUE)
CMH_TESTS        = CMH_TESTS[!is.na(CMH_TESTS$PVALUE),]

#--------------------------#
# 2) Compute FDR           #
#--------------------------#
CMH_TESTS$FDR = calculate_fdr(CMH_TESTS$PVALUE)

#--------------------------#
# 3) Save the result       #
#--------------------------#
saveRDS(CMH_TESTS, file="./data/tribolium_afc/CMH_TESTS_FDR.rds")
saveRDS(filter(CMH_TESTS, FDR<0.05), file="./data/tribolium_afc/CMH_TESTS_SIGNIFICANT.rds")

