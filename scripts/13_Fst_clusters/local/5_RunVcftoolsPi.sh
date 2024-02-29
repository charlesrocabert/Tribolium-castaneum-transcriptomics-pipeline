#!/bin/bash
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 5_RunVcftoolsPi.sh
# ------------------
# Run vcftools to compute the pi nucleotide diversity.
# (LOCAL SCRIPT)
#***************************************************************************

DATA_PATH="/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
ENVIRONMENTS=("CT" "HD")
GENERATIONS=("G1" "G21")
LINES=("L1" "L2" "L3" "L4" "L5" "L6" "Mx1" "Mx2")

# for env in ${ENVIRONMENTS[@]}; do
#   for gen in ${GENERATIONS[@]}; do
#     for line in ${LINES[@]}; do
#       echo "> Extracting VCF file for line "$env"_"$gen"_"$line
#       python $DATA_PATH/scripts/4_select_population/local/1_SelectPopulation.py -in-population ALL -in-suffix imputed -out-population $env\_$gen\_$line -out-suffix imputed -version Tcas3.30 -task VCF
#       mv $DATA_PATH/data/tribolium_snp/Tribolium_castaneum_$env\_$gen\_$line\_Tcas3.30_imputed.vcf $DATA_PATH/data/tribolium_snp/diversity/Tribolium_castaneum_$env\_$gen\_$line\_Tcas3.30_imputed.vcf
#     done
#   done
# done

# for env in ${ENVIRONMENTS[@]}; do
#   for gen in ${GENERATIONS[@]}; do
#     for line in ${LINES[@]}; do
#       echo "> Decompress GZ for line "$env"_"$gen"_"$line
#       GZVCF=$DATA_PATH/data/tribolium_snp/diversity/Tribolium_castaneum_$env\_$gen\_$line\_Tcas3.30_imputed.vcf.gz
#       if [ -s $GZVCF ]; then
#         gzip -d $GZVCF
#       else
#         echo "    No VCF file for this line"
#       fi
#     done
#   done
# done

# for env in ${ENVIRONMENTS[@]}; do
#   for gen in ${GENERATIONS[@]}; do
#     for line in ${LINES[@]}; do
#       echo "> Computing pi diverstity for line "$env"_"$gen"_"$line
#       VCF=$DATA_PATH/data/tribolium_snp/diversity/Tribolium_castaneum_$env\_$gen\_$line\_Tcas3.30_imputed.vcf
#       if [ -s $VCF ]; then
#         vcftools --vcf $VCF --window-pi  10000 --out $DATA_PATH/data/tribolium_diversity/pi/$env\_$gen\_$line
#       else
#         echo "    No VCF file for this line"
#       fi
#     done
#   done
# done

for env in ${ENVIRONMENTS[@]}; do
  for gen in ${GENERATIONS[@]}; do
    for line in ${LINES[@]}; do
      echo "> Compress to GZ for line "$env"_"$gen"_"$line
      VCF=$DATA_PATH/data/tribolium_snp/diversity/Tribolium_castaneum_$env\_$gen\_$line\_Tcas3.30_imputed.vcf
      if [ -s $VCF ]; then
        gzip $VCF
      else
        echo "    No VCF file for this line"
      fi
    done
  done
done
