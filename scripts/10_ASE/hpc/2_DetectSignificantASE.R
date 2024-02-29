#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_DetectSignificantASE.R
# ------------------------
# Detect significant ASE SNPs depending on population-level criteria.
# (HPC SCRIPT --> run wrapper)
#***************************************************************************

rm(list=ls())

### Load samples ###
load_samples <- function( population, version )
{
  filename = paste0("/scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_bam/samples_", population, "_", version, ".csv")
  samples = read.table(filename, h=T, sep=";")
  return(samples)
}

### Load ASE file ###
load_ASE_file <- function( population, version, suffix, sample )
{
  filename = paste0("/scratch/project_2003847/Tribolium_castaneum_ASE/tables/", population, "_", version, "_", suffix, "_", sample, ".table")
  ase_file = read.table(filename, h=T, sep="\t")
  return(ase_file)
}

### Calculate the FDR ###
calculate_fdr <- function( p )
{
  fdr = p.adjust(p, method="fdr")
  return(fdr)
}

### Run exact binomial tests per SNP and calculate FDR ###
run_exact_binomial_test <- function( ase_data )
{
  res = apply(X=ase_data,
              MARGIN = 1,
              FUN = function(t)
              {
                if (as.numeric(t[6]) > 0 & as.numeric(t[7]) > 0)
                {
                  ci = binom.test(x=as.numeric(t[7]), n=as.numeric(t[6])+as.numeric(t[7]), alternative="two.sided")
                  return(ci$p.value)
                }
                else
                {
                  return(NA)
                }
              })
  ase_data$p.value = res
  ase_data$fdr     = calculate_fdr(ase_data$p.value)
  return(ase_data)
}


##################
#      MAIN      #
##################

#-----------------------------------#
# 1) Read command line arguments    #
#-----------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
POPULATION = args[1]
VERSION    = args[2]
SUFFIX     = args[3]
PATH       = args[4]

setwd(PATH)

#-----------------------------------#
# 2) Load the list of samples       #
#-----------------------------------#
samples = load_samples(POPULATION, VERSION)

#-----------------------------------#
# 3) Load the complete list of SNPs #
#-----------------------------------#
snp_id = c()
for (sample in samples$sample)
{
  ase_data = load_ASE_file(POPULATION, VERSION, SUFFIX, sample)
  if (length(snp_id)==0)
  {
    snp_id = c(snp_id, ase_data$variantID)
  }
  else
  {
    snp_id = intersect(snp_id, ase_data$variantID)
  }
}
is_het_table              = data.frame(snp_id)
rownames(is_het_table)    = is_het_table$snp_id
is_signif_table           = data.frame(snp_id)
rownames(is_signif_table) = is_signif_table$snp_id

#-----------------------------------#
# 4) Compute ASE tests and FDR      #
#-----------------------------------#
for (sample in samples$sample)
{
  print(paste0(">>>> ", sample))
  ### 4.1) Load ASE data ###
  ase_data          = load_ASE_file(POPULATION, VERSION, SUFFIX, sample)
  ase_data          = run_exact_binomial_test(ase_data)
  ase_data$snp_id   = ase_data$variantID
  ase_data$isHet    = rep(0, dim(ase_data)[1])
  ase_data$isSignif = rep(0, dim(ase_data)[1])
  ### 4.2) Detect significant SNPs ###
  ase_data[ase_data$refCount>0 & ase_data$altCount>0,"isHet"]   = 1
  ase_data[!is.na(ase_data$fdr) & ase_data$fdr<0.05,"isSignif"] = 1
  rownames(ase_data) = ase_data$snp_id
  ### 4.3) Merge data ###
  is_het_table    = cbind(is_het_table, ase_data[is_het_table$snp_id,"isHet"])
  is_signif_table = cbind(is_signif_table, ase_data[is_signif_table$snp_id,"isSignif"])
}
names(is_het_table)[2:dim(is_het_table)[2]]       = samples$sample
names(is_signif_table)[2:dim(is_signif_table)[2]] = samples$sample

#-----------------------------------#
# 5) Extract significant SNPs based #
#    on population-level criteria   #
#-----------------------------------#
filename = paste0("/scratch/project_2003847/Tribolium_castaneum_ASE/", POPULATION, "_", VERSION, "_", SUFFIX, "_HET_SNPs")
write.table(is_het_table, file=paste0(filename,".txt"), quote=F, row.names=F, col.names=T)
saveRDS(is_het_table, file=paste0(filename,".rds"))

filename = paste0("/scratch/project_2003847/Tribolium_castaneum_ASE/", POPULATION, "_", VERSION, "_", SUFFIX, "_SIGNIF_SNPs")
write.table(is_signif_table, file=paste0(filename,".txt"), quote=F, row.names=F, col.names=T)
saveRDS(is_signif_table, file=paste0(filename,".rds"))

