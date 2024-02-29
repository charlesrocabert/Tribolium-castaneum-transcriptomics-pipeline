#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# plot_samples.R
# --------------
# Plot samples repartition in the experiment.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")
library("RColorBrewer")


##################
#      MAIN      #
##################

#------------------------------------#
# 1) Read command line arguments     #
#------------------------------------#
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3)
{
  stop("Please provide all command line arguments. Exit.", call.=FALSE)
}
PATH       = args[1]
POPULATION = args[2]
VERSION    = args[3]
#PATH       = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
#POPULATION = "ALL"
#VERSION    = "Tcas3.30"

setwd(PATH)

#------------------------------------#
# 2) Plot the repartition of samples #
#------------------------------------#
d = read.table(paste0("data/tribolium_bam/samples_",POPULATION,"_",VERSION,".csv"), h=T, sep=";")

ggplot(filter(d, source_env%in%c("CT","HD")), aes(source_env, fill=line)) + geom_bar(position='dodge') + facet_wrap(~ generation, ncol=1) + ggtitle("Repartition of individual samples in this experiment") + xlab("Condition") + ylab("Sampled individuals") + labs(fill="Line") + theme_minimal()

