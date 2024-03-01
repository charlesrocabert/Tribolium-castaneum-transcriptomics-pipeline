#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_PlotModuleExploration.R
# -------------------------
# Plot the result of the soft threshold exploration.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")


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
REPOSITORY_PATH = args[1]
POPULATION      = args[2]
VERSION         = args[3]
PHENOTYPE       = args[4]
THRESHOLD       = as.numeric(args[5])
setwd(REPOSITORY_PATH)

#--------------------------------#
# 2) Plot the result             #
#--------------------------------#
fitness_cor_data = read.table(paste0("./data/tribolium_modules/power_exploration_",POPULATION,"_",VERSION,"_",PHENOTYPE,".txt"), sep="\t", h=T)

p1 = ggplot(fitness_cor_data, aes(power, best_cor)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept=THRESHOLD, colour="Selected\nthreshold")) +
  xlab("Soft threshold (power)") +
  ylab("Corelation") +
  ggtitle("Module with best correlation to fitness") +
  theme_classic()
p2 = ggplot(fitness_cor_data, aes(power, log10(best_pval))) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept=THRESHOLD, colour="Selected\nthreshold")) +
  xlab("Soft threshold (power)") +
  ylab("P-value") +
  ggtitle("Associated p-value") +
  theme_classic()
p3 = ggplot(fitness_cor_data, aes(power, nb_modules)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept=THRESHOLD, colour="Selected\nthreshold")) +
  xlab("Number of modules") +
  ylab("P-value") +
  ggtitle("Number of modules found") +
  theme_classic()
p4 = ggplot(fitness_cor_data, aes(power, fit)) +
  geom_line() + geom_point() +
  geom_vline(aes(xintercept=THRESHOLD, colour="Selected\nthreshold")) +
  xlab("Number of modules") +
  ylab("Scale-free model fit") +
  ggtitle("Scale-free model fitting") +
  theme_classic()
legend = get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)) + scale_color_discrete(name="Legend"))
p = plot_grid(p1 + theme(legend.position="none"),
              p2 + theme(legend.position="none"),
              p3 + theme(legend.position="none"),
              p4 + theme(legend.position="none"),
              labels="AUTO")
p_legend = plot_grid(p, legend, rel_widths = c(3, .5))
ggsave(paste0("./data/tribolium_modules/module_exploration_",POPULATION,"_",VERSION,"_",PHENOTYPE,".pdf"), p_legend, width=9, height=7)

