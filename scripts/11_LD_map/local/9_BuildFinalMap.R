#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 9_BuildFinalMap.R
# -----------------
# Build the final genetic map and remove map ends at a given threshold.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")

plot_marey_map <- function( dataset, cut_thresholds, chromosome )
{
  p = ggplot(data=filter(dataset, chr==chromosome), aes(x=physical_pos, y=mean_pos)) +
    geom_point() +
    geom_hline(yintercept=cut_thresholds[chromosome][[1]][1], color="red") +
    geom_hline(yintercept=cut_thresholds[chromosome][[1]][2], color="red") +
    ggtitle(chromosome) +
    theme_minimal()
  return(p)
}

plot_marey_map_cut <- function( dataset, cut_thresholds, chromosome )
{
  # & mean_pos >= th1 & mean_pos <= th2 & physical_pos >= th3 & physical_pos <= th4
  th1 = cut_thresholds[chromosome][[1]][1]
  th2 = cut_thresholds[chromosome][[1]][2]
  th3 = cut_thresholds[chromosome][[1]][3]
  th4 = cut_thresholds[chromosome][[1]][4]
  p = ggplot() +
    geom_point(data=filter(dataset, chr==chromosome), aes(x=physical_pos, y=mean_pos), color="lightgrey") +
    geom_point(data=filter(dataset, chr==chromosome & mean_pos >= th1 & mean_pos <= th2 & physical_pos >= th3 & physical_pos <= th4), aes(x=physical_pos, y=mean_pos), color="blue") +
    geom_hline(yintercept=th1, color="red") +
    geom_hline(yintercept=th2, color="red") +
    geom_vline(xintercept=th3, color="red") +
    geom_vline(xintercept=th4, color="red") +
    ggtitle(chromosome) +
    theme_minimal()
  return(p)
}


##################
#      MAIN      #
##################

#---------------------------------------#
# 1) Read command line arguments        #
#---------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<6)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
SUFFIX          = args[4]
NUMBER_OF_LGS   = as.numeric(args[5])
THRESHOLD       = as.numeric(args[6])
setwd(REPOSITORY_PATH)

#---------------------------------------#
# 1) Define manual thresholds to cut    #
#    map ends                           #
#---------------------------------------#
cut_thresholds = list("ChLGX" = c(-1, 75, 0, 1.3e+7),
                      "ChLG2" = c(29, 87, 3e+6, 2e+7),
                      "ChLG3" = c(50, 145, 0, 4e+7),
                      "ChLG4" = c(42, 110, 0, 1.2e+7),
                      "ChLG5" = c(36, 101, 0, 2e+7),
                      "ChLG6" = c(37, 100, 0, 1.4e+7),
                      "ChLG7" = c(37, 100, 0, 2.2e+7),
                      "ChLG8" = c(48, 100, 0, 2.2e+7),
                      "ChLG9" = c(55, 120, 5e+6, 2.2e+7),
                      "ChLG10" = c(60, 125, 2.5e+6, 1.3e+7))

# plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG10")

#---------------------------------------#
# 2) Load each linkage group separately #
#    and build the genome-scale map     #
#---------------------------------------#
GENOME_SCALE_MAP = c()
for (lg in 1:NUMBER_OF_LGS)
{
  filename         = paste0("./data/tribolium_ld/",POPULATION,"_",VERSION,"_",SUFFIX,"_LG",lg,".txt")
  d                = read.table(filename, h=T, sep="\t")
  chr              = unique(d$chr)
  #dcut             = data.frame(d)
  dcut             = filter(d, mean_pos >= cut_thresholds[chr][[1]][1] & mean_pos <= cut_thresholds[chr][[1]][2] & physical_pos >= cut_thresholds[chr][[1]][3] & physical_pos <= cut_thresholds[chr][[1]][4])
  GENOME_SCALE_MAP = rbind(GENOME_SCALE_MAP, dcut)
}

#---------------------------------------#
# 3) Save tha dataset                   #
#---------------------------------------#
filename = paste0("./data/tribolium_ld/",POPULATION,"_",VERSION,"_",SUFFIX,"_ALL.txt")
write.table(GENOME_SCALE_MAP, file=filename, sep="\t", row.names=F, col.names =T, quote=F)

#---------------------------------------#
# 4) Plot the map                       #
#---------------------------------------#
# p1 = ggplot(GENOME_SCALE_MAP, aes(x=1:nrow(GENOME_SCALE_MAP), y=mean_pos)) +
#   geom_line() +
#   facet_wrap(~chr, scales="free") +
#   xlab("SNP index") +
#   ylab("Genetic position") +
#   theme_minimal()
#
# p2 = ggplot(GENOME_SCALE_MAP, aes(x=physical_pos, y=mean_pos)) +
#   geom_point() +
#   #geom_smooth(method="lm") +
#   facet_wrap(~chr, scales="free") +
#   xlab("Physical position") +
#   ylab("Genetic position") +
#   theme_minimal()
# p2
#
# pX = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLGX")
# p2 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG2")
# p3 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG3")
# p4 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG4")
# p5 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG5")
# p6 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG6")
# p7 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG7")
# p8 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG8")
# p9 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG9")
# p10 = plot_marey_map_cut(GENOME_SCALE_MAP, cut_thresholds, "ChLG10")
# plot_grid(pX, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=5)


