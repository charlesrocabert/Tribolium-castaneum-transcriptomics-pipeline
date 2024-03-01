#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_module_exploration_pipeline.sh
# --------------------------------
# Explore the soft power value to detect modules correlated to fitness.
# (LOCAL SCRIPT)
#***************************************************************************


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# The user must specify the path to the repository and GATK
REPOSITORY_PATH="/path/to/repository"
POPULATION="HD_G1"
VERSION="Tcas3.30"
COR_METHOD="pearson"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Run module exploration #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $REPOSITORY_PATH $POPULATION $VERSION $COR_METHOD EXPRESSION
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $REPOSITORY_PATH $POPULATION $VERSION $COR_METHOD PLASTICITY
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $REPOSITORY_PATH $POPULATION $VERSION $COR_METHOD NOISE

