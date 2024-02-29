#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_test_LODlimit_pipeline.sh
# ---------------------------
# Explore the LOD theshold for Lep-MAP3 SeparateChromosomes2
# (HPC SCRIPT --> run wrapper)
#***************************************************************************


for threshold in {5..20}
do
  sbatch python_lepMAP3_wrapper.sh SeparateChromosomes2.py -population CT_HD_G1 -version Tcas3.30 -lodlimit $threshold
done

