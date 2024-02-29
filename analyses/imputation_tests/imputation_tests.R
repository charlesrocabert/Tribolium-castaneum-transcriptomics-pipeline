#!/usr/bin/env Rscript

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# plot_imputation_tests.R
# -----------------------
# Plot the result of imputation tests.
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
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)<1)
# {
#   stop("Please provide all command line arguments. Exit.", call.=FALSE)
# }
# PATH = args[1]
PATH = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"

setwd(PATH)

#--------------------------------------#
# 2) Plot the probability distribution #
#    of missing genotypes              #
#--------------------------------------#
d = read.table("./data/tribolium_snp/imputation_tests/missing_genotypes_probability_distribution.csv", sep=";", h=F)
names(d) = "Nmissing"
d[,1] = d[,1]/822
ggplot(d, aes(Nmissing)) +
  geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
  xlab("% of missing genotypes per marker (1 - CR)") +
  ylab("Count") +
  ggtitle("Distribution of missing genotypes per marker") +
  theme_classic()

#--------------------------------------#
# 3) Plot the analysis                 #
#--------------------------------------#
d        = read.table("./data/tribolium_snp/imputation_tests/imputation_success_rate.csv", sep=";", h=T)
names(d) = c("success_rate", "positive", "total")
d        = d[1:100,]
N        = dim(d)[1]

ggplot(d, aes(x=success_rate)) +
  geom_histogram(fill="cornflowerblue", alpha=0.5, color="black") +
  xlab("") +
  ylab("Success rate") +
  ggtitle("Distribution of imputation\nsuccess rates over 100 repetitions") +
  theme_classic()

# ggplot(d, aes(x=total-positive)) +
#   geom_histogram() +
#   xlab("") +
#   ylab("Success rate") +
#   ggtitle("Distribution of success rates over 25 repetitions") +
#   theme_classic()

print(paste0("> Success rate mean = ", mean(d$success_rate), " (N=", N, ")"))
print(paste0("> Success rate se   = ", sd(d$success_rate)/N))

print(paste0("> Nb failures mean = ", mean(d$total-d$positive), " (N=", N, ")"))
print(paste0("> Nb failures se   = ", sd(d$total-d$positive)/N))

