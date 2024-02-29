#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 4_PlotAFCDistributions.R
# ------------------------
# Plot AFC distributions.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")

##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/")

#------------------------------------#
# 1) Load the AFC datasets  #
#------------------------------------#
AFC_CT = readRDS("./data/tribolium_afc/AFC_CT.rds")
AFC_HD = readRDS("./data/tribolium_afc/AFC_HD.rds")

ggplot() +
  geom_density(data=filter(AFC_CT, AFC!=0.0), aes(abs(AFC), color="CT")) +
  geom_density(data=filter(AFC_HD, AFC!=0.0), aes(abs(AFC), color="HD")) +
  scale_y_log10()

ggplot() +
  geom_density(data=AFC_CT, aes(abs(AFC), color="CT")) +
  geom_density(data=AFC_HD, aes(abs(AFC), color="HD")) +
  scale_y_log10()
