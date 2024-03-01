#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# read_counts_preparation_pipeline.sh
# -----------------------------------
# Read counts preparation pipeline.
# (LOCAL SCRIPT)
#***************************************************************************

VERSION="Tcas3.30"

### 1) Select populations from ALL read counts ###
python ../3_select_population/local/1_SelectPopulation.py -in-population ALL -in-suffix read_counts -out-population CT_HD_G1 -out-suffix read_counts -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectPopulation.py -in-population ALL -in-suffix read_counts -out-population CT_G1 -out-suffix read_counts -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectPopulation.py -in-population ALL -in-suffix read_counts -out-population HD_G1 -out-suffix read_counts -version $VERSION -task COUNTS

### 2) Transform CT-HD-G1 read counts ###
Rscript ./local/3_TransformReadCounts.R CT_HD_G1 read_counts CT_HD_G1 read_counts_READY $VERSION

### 3) Select sub-populations from CT/HD G1 transformed data ###
python ../3_select_population/local/1_SelectSubPopulation.py -in-population CT_HD_G1 -in-suffix read_counts_READY -out-population CT_G1 -out-suffix read_counts_transformed -version $VERSION -task COUNTS
python ../3_select_population/local/1_SelectSubPopulation.py -in-population CT_HD_G1 -in-suffix read_counts_READY -out-population HD_G1 -out-suffix read_counts_transformed -version $VERSION -task COUNTS

### 4) Detect low expressed transcripts ###
Rscript ./local/4_DetectLowExpressedTranscripts.R CT_G1 read_counts $VERSION
Rscript ./local/4_DetectLowExpressedTranscripts.R HD_G1 read_counts $VERSION

### 5) Standardize read counts ###
Rscript ./local/5_StandardizeReadCounts.R CT_G1 read_counts_transformed CT_G1 read_counts_READY $VERSION
Rscript ./local/5_StandardizeReadCounts.R HD_G1 read_counts_transformed HD_G1 read_counts_READY $VERSION

