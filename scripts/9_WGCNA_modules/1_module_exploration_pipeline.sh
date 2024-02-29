#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_module_exploration_pipeline.sh
# --------------------------------
# Explore the soft power value to detect modules correlated to fitness.
# (LOCAL SCRIPT)
#***************************************************************************


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
POPULATION="HD_G1"
VERSION="Tcas3.30"
COR_METHOD="pearson" # "pearson" / "bicor"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Run module exploration #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $POPULATION $VERSION $COR_METHOD EXPRESSION
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $POPULATION $VERSION $COR_METHOD PLASTICITY
Rscript ./local/1_ExploreModuleFitnessCorrelations.R $POPULATION $VERSION $COR_METHOD NOISE

