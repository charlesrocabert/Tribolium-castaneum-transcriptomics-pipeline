#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# LODlimit_exploration_merge_pipeline.sh
# --------------------------------------
# Merge map datasets from the LODlimit exploration
# (LOCAL SCRIPT)
#***************************************************************************

for val in {5..20}
do
  python ../MergeMapData.py -pcall CT_HD_G1_Tcas3.30_parent_call_filtered.txt -map LODlimit_exploration/CT_HD_G1_Tcas3.30_LOD$val\_map.txt -output LODlimit_exploration/LOD$val\_merged_map.txt
done


