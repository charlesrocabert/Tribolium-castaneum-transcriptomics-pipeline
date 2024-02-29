#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# variant_detection_pipeline.sh
# -----------------------------
# Variant detection pipeline to output ALL raw SNPs.
# (LOCAL SCRIPT)
#***************************************************************************

python ./local/4_SelectFilterAnnotateVariants.py -population ALL -version Tcas3.30 -filter-suffix raw_SNP

