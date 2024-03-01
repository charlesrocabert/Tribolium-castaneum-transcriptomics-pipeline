#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 8_MergeFinalMapData.py
# ----------------------
# Merge the final (ordered) map data (map + positions).
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
    parser.add_argument("--initial-map", "-initial-map", help="Initial map suffix")
    parser.add_argument("--ordered-map", "-ordered-map", help="Ordered map suffix")
    parser.add_argument("--parent-call", "-parent-call", help="Parent call suffix")
    parser.add_argument("--output", "-output", help="Output suffix")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Delete a file ###
def delete_file( filename ):
    os.system("rm "+filename)
    assert_deletion(filename)

### Load the initial map and its physical positions ###
def load_initial_map( population, version, initial_map_suffix, parent_call_suffix ):
    map_file     = "./data/tribolium_ld/"+population+"_"+version+"_"+initial_map_suffix+".txt"
    pcall_file   = "./data/tribolium_ld/"+population+"_"+version+"_"+parent_call_suffix+".txt"
    chrom_length = {"ChLG10": 11385945,
                    "ChLG2": 20218313,
                    "ChLG3": 38791426,
                    "ChLG4": 13894384,
                    "ChLG5": 19135674,
                    "ChLG6": 13176117,
                    "ChLG7": 20532118,
                    "ChLG8": 18021249,
                    "ChLG9": 21459566,
                    "ChLGX": 10876780}
    #----------------------------------------#
    # 1) Open files                          #
    #----------------------------------------#
    map   = open(map_file, "r")
    pcall = open(pcall_file, "r")
    #----------------------------------------#
    # 2) Manage headers                      #
    #----------------------------------------#
    for i in range(1):
        map.readline()
    for i in range(7):
        pcall.readline()
    #----------------------------------------#
    # 3) Load the dataset                    #
    #----------------------------------------#
    markers = {}
    l1      = map.readline()
    l2      = pcall.readline()
    count   = 0
    while l1:
        assert l2
        if count % 10000 == 0:
            print(">> "+str(count)+" markers parsed")
        count  += 1
        l1      = l1.strip("\n").split("\t")
        l2      = l2.strip("\n").split("\t")
        lg      = l1[0]
        chr     = l2[0]
        pos     = l2[1]
        chr_len = chrom_length[chr]
        assert count not in markers.keys()
        markers[count] = {"lg":lg, "chr":chr, "pos":pos, "chr_len":str(chr_len)}
        l1 = map.readline()
        l2 = pcall.readline()
    assert not l2
    map.close()
    pcall.close()
    print(">> "+str(count)+" markers parsed in total.")
    #----------------------------------------#
    # 4) Return the complete list of markers #
    #----------------------------------------#
    return markers

### Load the ordered map and its genomic positions ###
def load_ordered_map( population, version, ordered_map_suffix ):
    #-------------------------------#
    # 1) Open files                 #
    #-------------------------------#
    map_file = "./data/tribolium_ld/"+population+"_"+version+"_"+ordered_map_suffix+".txt"
    map      = open(map_file, "r")
    #-------------------------------#
    # 2) Load the dataset           #
    #-------------------------------#
    genomic_positions = {}
    l                 = map.readline()
    count             = 0
    LG                = 0
    likelihood        = 0.0
    while l:
        if l.startswith("#"):
            if l.startswith("#*** LG = "):
                l          = l.strip("\n#* ").split(" ")
                lg         = int(l[2])
                likelihood = float(l[5])
        else:
            if count % 10000 == 0:
                print(">> "+str(count)+" markers parsed")
            count     += 1
            l          = l.strip("\n").split("\t")
            number     = int(l[0])
            male_pos   = l[1]
            female_pos = l[2]
            assert number not in genomic_positions.keys()
            genomic_positions[number] = {"male_pos":male_pos, "female_pos":female_pos, "lg":str(lg), "likelihood":str(likelihood)}
        l = map.readline()
    map.close()
    print(">> "+str(count)+" markers parsed in total.")
    #-------------------------------#
    # 3) Return the list of markers #
    #-------------------------------#
    return genomic_positions

### Merge physical and genomic informations ###
def merge_maps( population, version, output_suffix, markers, genomic_positions ):
    f1 = open("./data/tribolium_ld/"+population+"_"+version+"_"+output_suffix+".txt", "w")
    #f2 = open("./data/tribolium_ld/Tribolium_castaneum_"+version+"_plink_map.txt", "w")
    f1.write("lg\tchr\tphysical_pos\tchr_len\tmale_pos\tfemale_pos\tmean_pos\tlikelihood\n")
    for item in genomic_positions.items():
        assert item[0] in markers.keys(), item[0]
        assert markers[item[0]]["lg"] == genomic_positions[item[0]]["lg"]
        marker = markers[item[0]]
        ######################
        line1  = ""
        line1 += marker["lg"]+"\t"
        line1 += marker["chr"]+"\t"
        line1 += marker["pos"]+"\t"
        line1 += marker["chr_len"]+"\t"
        line1 += item[1]["male_pos"]+"\t"
        line1 += item[1]["female_pos"]+"\t"
        line1 += str((float(item[1]["male_pos"])+float(item[1]["female_pos"]))/2)+"\t"
        line1 += item[1]["likelihood"]+"\n"
        ######################
        # line2  = ""
        # line2 += marker["chr"]+"\t"
        # line2 += marker["chr"]+"_"+marker["pos"]+"\t"
        # line2 += item[1]["female_pos"]+"\t"
        # line2 += marker["pos"]+"\n"
        ######################
        f1.write(line1)
        #f2.write(line2)
    f1.close()
    #f2.close()


##################
#      MAIN      #
##################

if __name__ == '__main__':
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments                    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(config["repository_path"])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load the initial map and its physical positions #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the initial map and its physical positions")
    markers = load_initial_map(config["population"], config["version"], config["initial_map"], config["parent_call"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Load the ordered map and its genomic positions  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load the ordered map and its genomic positions")
    genomic_positions = load_ordered_map(config["population"], config["version"], config["ordered_map"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Merge maps in a single file                     #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Merge maps in a single file")
    merge_maps(config["population"], config["version"], config["output"], markers, genomic_positions)

    print(">> Done.")

