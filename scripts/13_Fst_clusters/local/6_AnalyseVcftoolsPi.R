#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 6_AnalyseVcftoolsPi.R
# ---------------------
# Analyse vcftools pi diversity measures and generate figures.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("reshape2")
library("corrplot")
library("ggcorrplot")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

# ENVIRONMENTS = c("CT", "HD")
# GENERATIONS  = c("G1", "G21")
# LINES        = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2")

ENVIRONMENTS = c("HD")
GENERATIONS  = c("G1", "G21")
LINES        = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2")

D = data.frame()
for (env in ENVIRONMENTS)
{
  for (gen in GENERATIONS)
  {
    for (line in LINES)
    {
      filename = paste0("data/tribolium_diversity/pi/",env,"_",gen,"_",line,".windowed.pi")
      print(filename)
      if (file.exists(filename))
      {
        d    = read.table(filename, h=T, sep="\t")
        name = paste0(line,"-",env,"-",gen)
        res  = data.frame(rep(name, dim(d)[1]), d$PI)
        names(res) = c("line", "pi")
        D    = rbind(D, res)
      }
    }
  }
}

ggplot(D, aes(line, log10(pi))) + geom_boxplot()




