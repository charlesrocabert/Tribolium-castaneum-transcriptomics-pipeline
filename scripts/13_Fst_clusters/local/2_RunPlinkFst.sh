#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 2_RunPlinkFst.sh
# ----------------
# Run PLINK 2.0 to compute pairwise population Fst.
# (LOCAL SCRIPT)
#***************************************************************************

/Users/charlesrocabert/plink2/plink2 -vcf ../../../data/tribolium_snp/Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf -fst pop method="hudson" -pheno ../../../data/tribolium_diversity/fst/CT_HD_Tcas3.30_clusters.txt --allow-extra-chr --out ../../../data/tribolium_diversity/fst/CT_HD_Tcas3.30_HUDSON

/Users/charlesrocabert/plink2/plink2 -vcf ../../../data/tribolium_snp/Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf -fst pop method="wc" -pheno ../../../data/tribolium_diversity/fst/CT_HD_Tcas3.30_clusters.txt --allow-extra-chr --out ../../../data/tribolium_diversity/fst/CT_HD_Tcas3.30_WC

