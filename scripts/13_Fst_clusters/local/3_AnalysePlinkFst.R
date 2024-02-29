#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_AnalysePlinkFst.R
# -------------------
# Analyse PLINK 2.0 Fst results and generate figures.
# (LOCAL SCRIPT)
#***************************************************************************

library("tidyverse")
library("reshape2")
library("corrplot")
library("ggcorrplot")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")

### Rebuild the dataset ###
rebuild_dataset <- function( fst_data )
{
  ### Rebuild ###
  line1       = c()
  env1        = c()
  generation1 = c()
  pop1        = c()
  line2       = c()
  env2        = c()
  generation2 = c()
  pop2        = c()
  fst         = c()
  for(i in 1:dim(fst_data)[1])
  {
    elmts1 = strsplit(d[i,1], "-")
    elmts2 = strsplit(d[i,2], "-")
    # ### Right column ###
    # line1       = c(line1, elmts1[[1]][1], elmts2[[1]][1])
    # env1        = c(env1, elmts1[[1]][2], elmts2[[1]][2])
    # generation1 = c(generation1, elmts1[[1]][3], elmts2[[1]][3])
    # pop1        = c(pop1, d[i,1], d[i,2])
    # ### Left column ###
    # line2       = c(line2, elmts2[[1]][1], elmts1[[1]][1])
    # env2        = c(env2, elmts2[[1]][2], elmts1[[1]][2])
    # generation2 = c(generation2, elmts2[[1]][3], elmts1[[1]][3])
    # pop2        = c(pop2, d[i,2], d[i,1])
    # ### Fst ###
    # fst = c(fst, d[i,3], d[i,3])
    ### Right column ###
    line1       = c(line1, elmts1[[1]][1])
    env1        = c(env1, elmts1[[1]][2])
    generation1 = c(generation1, elmts1[[1]][3])
    pop1        = c(pop1, d[i,1])
    ### Left column ###
    line2       = c(line2, elmts2[[1]][1])
    env2        = c(env2, elmts2[[1]][2])
    generation2 = c(generation2, elmts2[[1]][3])
    pop2        = c(pop2, d[i,2])
    ### Fst ###
    fst = c(fst, d[i,3])
  }
  newd = data.frame(cbind(pop1, line1, env1, generation1, line2, pop2, env2, generation2, fst))
  newd$generation1 = as.numeric(newd$generation1)
  newd$generation2 = as.numeric(newd$generation2)
  newd$fst         = as.numeric(newd$fst)
  # newd = data.frame(cbind(pop1, line1, env1, generation1, fst))
  # newd$generation1 = as.numeric(newd$generation1)
  # newd$fst         = as.numeric(newd$fst)
  return(newd)
}

### Extract Fst values for a given line ###
extract_fst <- function( data, env, gen, line )
{
  data = filter(data, env1==env & env2==env & generation1==gen & generation2==gen)
  data = filter(data, (env1==env & generation1==gen & line1==line) | (env2==env & generation2==gen & line2==line))
  if (dim(data)[1] > 0)
  {
    #res = data.frame(rep(paste0(env,"_",gen,"_",line), dim(data)[1]), data$fst)
    res = data.frame(rep(paste0(line), dim(data)[1]), data$fst)
    names(res) = c("line", "fst")
    return(res)
  } else {
    return(F)
  }
}


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

### Load FST ###
d = read.table("data/tribolium_diversity/fst/CT_HD_Tcas3.30_WC.fst.summary", h=T, sep="\t")

### Rebuild dataset ###
d = rebuild_dataset(d)

### Extract the list of Fst for each line ###
# ENVIRONMENTS = c("CT", "HD")
# GENERATIONS  = c(1, 21)
# LINES        = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2")
ENVIRONMENTS = c("HD")
GENERATIONS  = c(21)
LINES        = c("L1", "L2", "L3", "L4", "L5", "L6", "Mx1", "Mx2")

D = data.frame()
for (env in ENVIRONMENTS)
{
  for (gen in GENERATIONS)
  {
    for (line in LINES)
    {
      res = extract_fst(d, env, gen, line)
      if (is.data.frame(res))
      {
        D = rbind(D, res)
      }
    }
  }
}

# my_comp = list(c("HD_21_L2", "HD_21_L1"),
#                c("HD_21_L2", "HD_21_L3"),
#                c("HD_21_L2", "HD_21_L5"),
#                c("HD_21_L2", "HD_21_L6"),
#                c("HD_21_L2", "HD_21_Mx1"),
#                c("HD_21_L2", "HD_21_Mx2"))
my_comp = list(c("L2", "L1"),
               c("L2", "L3"),
               c("L2", "L5"),
               c("L2", "L6"),
               c("L2", "Mx1"),
               c("L2", "Mx2"))

ggplot(D, aes(line, fst, fill=line)) +
  geom_boxplot() +
  scale_fill_brewer(palette="BrBG") +
  stat_compare_means(comparisons=my_comp, method="wilcox.test", vjust=-0.5) +
  ggtitle("") +
  xlab("Line") +
  ylab("Pairwise Fst") +
  theme_classic()

