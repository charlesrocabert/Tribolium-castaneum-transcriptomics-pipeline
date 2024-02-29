#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_CheckBamReadgroups.py
# -----------------------
# Check the absence of the readgroup in every Tcas3.30 (2016) and Tcas5.2
# (2017) BAM files.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import subprocess

### Load the bam map ###
def load_bam_map( filename ):
    samples   = []
    file      = open(filename, "r")
    csvreader = csv.reader(file, delimiter=";")
    header    = next(csvreader)
    samples.append(header)
    for row in csvreader:
        sample = {}
        for i in range(len(header)):
            sample[header[i]] = row[i]
        samples.append(sample)
    file.close()
    return samples

### Check the existence of a read group ###
def check_readgroup( samples ):
    for sample in samples[1:]:
        bam_path = sample["bam_path"]
        proc     = subprocess.Popen(["samtools view -H "+bam_path+" | grep ^\@RG"], stdout=subprocess.PIPE, shell=True)
        output   = proc.stdout.read().decode('utf8')
        if len(output.strip("\n")) > 0:
            print("> sample "+sample["sample"]+" "+sample["genome_ref"]+" "+sample["annotation"])
            print("  >>> Read group is present.")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 3_CheckBamReadgroups.py")
    print("# -----------------------")
    print("# Check the absence of the readgroup in every Tcas3.30 (2016) and Tcas5.2")
    print("# (2017) BAM files.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    DD_PATH = "/Volumes/TRIBOLIUM"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Check Tcas3.30 BAM read groups #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print("> Evaluate Tcas3.30 (2016) samples")
    tcas3_30_samples = load_bam_map("data/tribolium_bam/bam_map_ALL_Tcas3.30.csv")
    check_readgroup(tcas3_30_samples)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Check Tcas5.2 BAM read groups  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print("> Evaluate Tcas5.2 (2017) samples")
    tcas5_2_samples = load_bam_map("data/tribolium_bam/bam_map_ALL_Tcas5.2.csv")
    check_readgroup(tcas5_2_samples)

