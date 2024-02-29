#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_eQTLs_data_preparation_pipeline.sh
# ------------------------------------
# eQTLs data preparation pipeline (before eQTLs mapping).
# (LOCAL SCRIPT)
#***************************************************************************


cd /Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1) Define main parameters  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
VERSION="Tcas3.30"
IMPUTATION="_imputed"
PLINK2="/Users/charlesrocabert/plink2/plink2"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2) Build gemma input files #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Generate .bed and .bim files"
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_CT_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/CT_G1_$VERSION$IMPUTATION\_EXPRESSION
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_CT_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/CT_G1_$VERSION$IMPUTATION\_FITNESS
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_EXPRESSION
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_PLASTICITY
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_NOISE
$PLINK2 --vcf ./data/tribolium_snp/Tribolium_castaneum_HD_G1_$VERSION\_eQTL$IMPUTATION.vcf --allow-extra-chr --make-bed --out ./data/tribolium_eqtl/gemma/HD_G1_$VERSION$IMPUTATION\_FITNESS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3) Edit Fam files          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo ">> Edit .fam files"
Rscript ./scripts/10_eQTLs/local/EditFam.R CT_G1 $VERSION $IMPUTATION EXPRESSION
Rscript ./scripts/10_eQTLs/local/EditFam.R CT_G1 $VERSION $IMPUTATION FITNESS
Rscript ./scripts/10_eQTLs/local/EditFam.R HD_G1 $VERSION $IMPUTATION EXPRESSION
Rscript ./scripts/10_eQTLs/local/EditFam.R HD_G1 $VERSION $IMPUTATION PLASTICITY
Rscript ./scripts/10_eQTLs/local/EditFam.R HD_G1 $VERSION $IMPUTATION NOISE
Rscript ./scripts/10_eQTLs/local/EditFam.R HD_G1 $VERSION $IMPUTATION FITNESS

