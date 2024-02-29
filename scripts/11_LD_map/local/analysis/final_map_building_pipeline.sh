#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# final_map_building_pipeline.sh
# ------------------------------
# Build the final genetic map datasets.
# (LOCAL SCRIPT)
#***************************************************************************

for lg in {1..10}
do
  python ../8_MergeFinalMapData.py -population CT_HD_G1 -version Tcas3.30 -initial-map map_cleaned -ordered-map map_ordered_LG$lg -parent-call parent_call_cleaned -output final_map_LG$lg
done

