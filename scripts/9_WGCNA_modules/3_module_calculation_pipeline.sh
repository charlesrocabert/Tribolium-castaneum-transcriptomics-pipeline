#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_module_calculation_pipeline.sh
# --------------------------------
# Build modules and save modules correlated to fitness.
# (LOCAL SCRIPT)
#***************************************************************************


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
POPULATION="HD_G1"
VERSION="Tcas3.30"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Run module exploration #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rscript ./local/3_ComputeModules.R $POPULATION $VERSION EXPRESSION 8
Rscript ./local/3_ComputeModules.R $POPULATION $VERSION PLASTICITY 7
Rscript ./local/3_ComputeModules.R $POPULATION $VERSION NOISE 6

