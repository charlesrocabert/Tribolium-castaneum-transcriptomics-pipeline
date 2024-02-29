#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# phenotype_preparation_pipeline.sh
# ---------------------------------
# Phenotypes preparation pipeline.
# (LOCAL SCRIPT)
#***************************************************************************

VERSION="Tcas3.30"

Rscript ./local/1_CopyExpressionFiles.sh
Rscript ./local/2_ComputePlasticityResponse.R
Rscript ./local/3_ComputePhenotypicNoise.R
Rscript ./local/4_ComputeRelativeFitness.R

