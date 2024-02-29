#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_ordermarkers2_pipeline.sh
# ---------------------------
# Run Lep-MAP3 OrderMarkers2, chromosome per chromosome.
# (HPC SCRIPT --> run wrapper)
#***************************************************************************


for chr in {1..10}
do
  sbatch python_lepMAP3_wrapper.sh OrderMarkers2.py -population CT_HD_G1 -version Tcas3.30 -in-pcall parent_call_cleaned -in-map map_cleaned -out-map map_ordered_LG$chr -chr $chr
done

