#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_plot_module_exploration_pipeline.sh
# -------------------------------------
# Plot the soft power value exploration.
# (LOCAL SCRIPT)
#***************************************************************************


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
POPULATION="HD_G1"
VERSION="Tcas3.30"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Run module exploration #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rscript ./local/2_PlotModuleExploration.R $REPOSITORY_PATH $POPULATION $VERSION EXPRESSION 8
Rscript ./local/2_PlotModuleExploration.R $REPOSITORY_PATH $POPULATION $VERSION PLASTICITY 7
Rscript ./local/2_PlotModuleExploration.R $REPOSITORY_PATH $POPULATION $VERSION NOISE 6

