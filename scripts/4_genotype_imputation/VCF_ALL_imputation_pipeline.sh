#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# VCF_ALL_imputation_pipeline.sh
# ------------------------------
# Genotypes imputation pipeline (ALL population).
# (LOCAL SCRIPT)
#***************************************************************************

python ./local/3_ImputeGenotypes.py -population ALL -version Tcas3.30 -suffix Beagle

