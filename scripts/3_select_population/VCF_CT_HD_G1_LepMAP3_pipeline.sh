#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# VCF_CT_HD_G1_LepMAP3_pipeline.sh
# --------------------------------
# SNP pipeline for Lep-MAP3 (CT/HD G1 population).
# (LOCAL SCRIPT)
#***************************************************************************

python ./local/1_SelectPopulation.py -in-population ALL -in-suffix raw_SNP -out-population CT_HD_G1 -out-suffix raw_SNP -version Tcas3.30 -task VCF
python ./local/2_FilterGenotypes.py -population CT_HD_G1 -version Tcas3.30 -suffix LepMAP3
python ./local/3_VariantsToTable.py -population CT_HD_G1 -version Tcas3.30 -suffix LepMAP3 -annotation

