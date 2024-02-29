#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 4_RelocateBamFiles.py
# ---------------------
# Relocate the BAM files in a specific folder for each version.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import time
import subprocess

### Load the BAM map ###
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

### Write the BAM map in a specified file ###
def write_bam_map( samples, filename ):
    f      = open(filename, "w")
    header = samples[0]
    f.write(";".join(header)+"\n")
    for sample in samples[1:]:
        line = []
        for elmt in header:
            line.append(sample[elmt])
        f.write(";".join(line)+"\n")
    f.close()

### Initialize the relocation folder ###
def init_folder( path, version ):
    # Confirm with the user
    print("> You are about to erase the folder "+path+"/"+version)
    user_input = input('Confirm? [Y/N] ')
    # Input validation
    if user_input.lower() not in ('y', 'yes'):
        print("> Exit.")
        sys.exit()
    # Delete data
    if os.path.exists(PATH+"/"+FOLDER):
        os.system("rm -rf "+PATH+"/"+FOLDER)
    os.system("mkdir "+PATH+"/"+FOLDER)

### Relocate the files ###
def relocate_bam( samples, path, version ):
    relocated_samples = []
    relocated_samples.append(samples[0])
    for sample in samples[1:]:
        new_bam_path = path+"/"+version+"/"+sample["sample"]+".bam"
        new_bai_path = path+"/"+version+"/"+sample["sample"]+".bai"
        if not os.path.isfile(new_bam_path):
            print("> Copy BAM "+sample["sample"])
            os.system("cp "+sample["bam_path"]+" "+new_bam_path)
        sample["bam_path"] = new_bam_path
        sample["bai_path"] = new_bai_path
        relocated_samples.append(sample)
    write_bam_map(relocated_samples, "data/tribolium_bam/bam_map_ALL_"+version+"_relocated.csv")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 4_RelocateBamFiles.py")
    print("# ---------------------")
    print("# Relocate the BAM files in a specific folder for each version.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    DD_PATH = "/Volumes/TRIBOLIUM/TRIBOLIUM-BAM"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Relocate Tcas3.30 BAM files #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    version = "Tcas3.30"
    samples = load_bam_map("data/tribolium_bam/bam_map_ALL_"+version+".csv")
    init_folder(DD_PATH, V1)
    relocate_bam(samples, DD_PATH, V1)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Relocate Tcas5.2 BAM files  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    version = "Tcas5.2"
    samples = load_bam_map("data/tribolium_bam/bam_map_ALL_"+version+".csv")
    init_folder(DD_PATH, V2)
    relocate_bam(samples, DD_PATH, V2)

