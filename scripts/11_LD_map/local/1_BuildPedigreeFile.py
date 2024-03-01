#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 1_BuildPedigreeFile.py
# ----------------------
# Build a pedigree file from a sample list file.
# (LOCAL SCRIPT)
#***************************************************************************

import os
import sys
import csv
import time
import argparse
import subprocess

### Parse command line arguments ###
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-path", "-repository-path", help="Repository path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Parse the samples file ###
def parse_samples_file( population, version ):
    families = {}
    f        = open("./data/tribolium_bam/samples_"+population+"_"+version+".csv", "r")
    l        = f.readline()
    l        = f.readline()
    while l:
        l = l.strip("\n").split(";")
        progeny = l[0]
        family  = l[2]
        mother  = l[2]
        father  = l[5]
        if family not in families.keys():
            families[family] = {"family":family, "mother":mother, "father":father, "progeny":[progeny]}
        else:
            assert mother == families[family]["mother"]
            assert father == families[family]["father"]
            families[family]["progeny"].append(progeny)
        l = f.readline()
    f.close()
    return families

### Parse the GFF file ###
def build_pedigree_file( families, population, version ):
    f     = open("./data/tribolium_pedigree/pedigree_"+population+"_"+version+".txt", "w")
    ### Build line 1 ###
    line1 = "CHR\tPOS"
    line2 = "CHR\tPOS"
    line3 = "CHR\tPOS"
    line4 = "CHR\tPOS"
    line5 = "CHR\tPOS"
    line6 = "CHR\tPOS"
    for item in families.items():
        line1 += "\t"+item[1]["family"]+"\t"+item[1]["family"]
        line2 += "\t"+item[1]["mother"]+"\t"+item[1]["father"]
        line3 += "\t0\t0"
        line4 += "\t0\t0"
        line5 += "\t2\t1"
        line6 += "\t0\t0"
        for progeny in item[1]["progeny"]:
            line1 += "\t"+item[1]["family"]
            line2 += "\t"+progeny
            line3 += "\t"+item[1]["father"]
            line4 += "\t"+item[1]["mother"]
            line5 += "\t2"
            line6 += "\t0"
    f.write(line1+"\n")
    f.write(line2+"\n")
    f.write(line3+"\n")
    f.write(line4+"\n")
    f.write(line5+"\n")
    f.write(line6+"\n")
    f.close()


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Parse the samples file       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse the samples file")
    families = parse_samples_file(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Build the pedigree file      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Build the pedigree file")
    build_pedigree_file(families, config["population"], config["version"])

    print(">> Done.")

