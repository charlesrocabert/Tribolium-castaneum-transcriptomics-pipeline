#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# plot_ld_map.R
# -------------
# Plot a genetic map.
# (LOCAL SCRIPT)
#***************************************************************************

rm(list=ls())

library("tidyverse")
library("cowplot")


##################
#      MAIN      #
##################

setwd("/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation")

d = read.table("./data/tribolium_ld/CT_HD_G1_Tcas3.30_map_joined.txt", h=T, sep="\t")

head(d)



# d = read.table("./data/tribolium_ld/CT_HD_G1_Tcas3.30_final_map_ALL_0.1.txt", h=T, sep="\t")
#
# par(mfrow=c(2,5))
# for(CHR in unique(d$chr))
# {
#   print(CHR)
#   dl = filter(d, chr==CHR)
#   HP = sort(table(dl$mean_pos), decreasing=T)
#   #HP = HP[HP>=20]
#   barplot(HP, main=CHR)
#   legend("topright", legend=c(paste0(sum(HP>=20)," haploblocks with N >= 20")), bty="n")
#   abline(h=20, col="red")
# }
