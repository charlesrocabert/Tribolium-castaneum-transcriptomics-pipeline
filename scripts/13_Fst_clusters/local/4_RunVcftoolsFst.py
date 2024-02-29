#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_RunVcftoolsFst.py
# -------------------
# Run SNP-level FSTs for each line/environment/generation.
#***************************************************************************

import os
import sys
import csv

### Load the list of samples corresponding to a given sub-population ###
def extract_samples_list( environment, generation, line, pop_name ):
    samples_list = []
    filename     = "./data/tribolium_bam/samples_"+environment+"_"+generation+"_"+line+"_Tcas3.30.csv"
    if os.path.isfile(filename):
        file      = open(filename, "r")
        csvreader = csv.reader(file, delimiter=";")
        header    = next(csvreader)
        for row in csvreader:
            samples_list.append(row[0])
        file.close()
        if len(samples_list) > 0:
            file = open("./data/tribolium_diversity/fst/"+pop_name+".txt", "w")
            for sample in samples_list:
                file.write(sample+"\n")
            file.close()
            return True
        return False
    return False

### Run the calculation of SNP-level Fst with vcftools ###
def run_vcftools_fst( vcf_file, input, output ):
    cmdline = []
    cmdline.append("vcftools")
    cmdline.append("--vcf "+vcf_file)
    cmdline.append("--weir-fst-pop "+input)
    cmdline.append("--out ./data/tribolium_diversity/fst/"+output)
    os.system(" ".join(cmdline))


##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 3_RunVcftoolsFst.py")
    print("# -------------------")
    print("# Run SNP-level FSTs for each line/environment/generation.")
    print("#***************************************************************************")
    print("")

    WD_PATH = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    os.chdir(WD_PATH)

    ENVIRONMENTS = ["CT", "HD"]
    LINES        = ["L1", "L2", "L3", "L5", "L6", "Mx1", "Mx2"]
    GENERATIONS  = ["G1", "G21"]
    VCF         = "./data/tribolium_snp/Tribolium_castaneum_ALL_Tcas3.30_imputed.vcf"
    for env in ENVIRONMENTS:
        for gen in GENERATIONS:
            for line in LINES:
                print("> Compute SNP-level FST for environment "+env+", generation "+gen+" and line "+line)
                output_name = env+"_"+gen+"_"+line
                success = extract_samples_list(env, gen, line, output_name)
                if success:
                    run_vcftools_fst(VCF, output_name, output_name)

















