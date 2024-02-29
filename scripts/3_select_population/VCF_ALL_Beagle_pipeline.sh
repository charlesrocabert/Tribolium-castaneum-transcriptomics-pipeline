#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# VCF_ALL_Beagle_pipeline.sh
# --------------------------
# SNP pipeline for Beagle (ALL population).
# (LOCAL SCRIPT)
#***************************************************************************

python ./local/2_FilterGenotypes.py -population ALL -version Tcas3.30 -suffix Beagle
python ./local/3_VariantsToTable.py -population ALL -version Tcas3.30 -suffix Beagle -annotation

