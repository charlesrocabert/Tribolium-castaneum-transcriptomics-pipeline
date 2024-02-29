#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# VCF_HD_G1_eQTL_imputed_pipeline.sh
# ----------------------------------
# SNP pipeline for eQTL (HD G1 imputed population).
# (LOCAL SCRIPT)
#***************************************************************************

python ./local/1_SelectPopulation.py -in-population ALL -in-suffix imputed -out-population HD_G1 -out-suffix imputed -version Tcas3.30 -task VCF
python ./local/2_FilterGenotypes.py -population HD_G1 -version Tcas3.30 -imputed -suffix eQTL_imputed
python ./local/3_VariantsToTable.py -population HD_G1 -version Tcas3.30 -suffix eQTL_imputed -annotation

