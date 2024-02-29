#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 6_merge_eQTLs_pipeline.sh
# -------------------------
# Merge eQTL associations.
# (LOCAL SCRIPT)
#***************************************************************************


cd /Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Compute the correlation between phenotype and fitness #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Compute the correlation between phenotype and fitness"
Rscript ./scripts/10_eQTLs/local/ComputeCorrelationToFitness.R CT_G1 Tcas3.30 EXPRESSION
Rscript ./scripts/10_eQTLs/local/ComputeCorrelationToFitness.R HD_G1 Tcas3.30 EXPRESSION
Rscript ./scripts/10_eQTLs/local/ComputeCorrelationToFitness.R HD_G1 Tcas3.30 PLASTICITY
Rscript ./scripts/10_eQTLs/local/ComputeCorrelationToFitness.R HD_G1 Tcas3.30 NOISE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Extract the list of gene positions                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Extract the list of gene positions"
python ./scripts/10_eQTLs/local/ExtractGenePos.py -version Tcas3.30

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3) Merge eQTLs data with annotations                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Merge eQTLs data with annotations and fitness (imputed data)"
Rscript ./scripts/10_eQTLs/local/MergeEQTLsDatasets.R CT_G1 Tcas3.30 _imputed EXPRESSION CT
Rscript ./scripts/10_eQTLs/local/MergeEQTLsDatasets.R HD_G1 Tcas3.30 _imputed EXPRESSION HD
Rscript ./scripts/10_eQTLs/local/MergeEQTLsDatasets.R HD_G1 Tcas3.30 _imputed PLASTICITY HD
Rscript ./scripts/10_eQTLs/local/MergeEQTLsDatasets.R HD_G1 Tcas3.30 _imputed NOISE HD

