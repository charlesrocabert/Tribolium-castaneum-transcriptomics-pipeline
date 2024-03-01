#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# MergeMapData.py
# ---------------
# Merge map data (map + positions).
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
    parser.add_argument("--pcall", "-pcall", help="Parent call filename")
    parser.add_argument("--map", "-map", help="Map filename")
    parser.add_argument("--output", "-output", help="Output filename")
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

### Merge map data with position information ###
def merge_map_data( pcall_filename, map_filename, output_filename ):
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
    #-------------------#
    # 1) Open files     #
    #-------------------#
    pcall  = open("./data/tribolium_ld/"+pcall_filename, "r")
    map    = open("./data/tribolium_ld/"+map_filename, "r")
    output = open("./data/tribolium_ld/"+output_filename, "w")
    output.write("LG\tCHR\tPOS\tCHR_LEN\n")
    output.flush()
    #-------------------#
    # 2) Manage headers #
    #-------------------#
    for i in range(1):
        map.readline()
    for i in range(7):
        pcall.readline()
    #-------------------#
    # 3) Merge datasets #
    #-------------------#
    l1    = map.readline()
    l2    = pcall.readline()
    count = 0
    while l1:
        assert l2
        if count % 10000 == 0:
            print(">> "+str(count)+" markers parsed")
        count += 1
        l1  = l1.strip("\n").split("\t")
        l2  = l2.strip("\n").split("\t")
        lg  = l1[0]
        chr = l2[0]
        pos = l2[1]
        output.write(lg+"\t"+chr+"\t"+pos+"\t"+str(chrom_length[chr])+"\n")
        output.flush()
        l1 = map.readline()
        l2 = pcall.readline()
    assert not l2
    map.close()
    pcall.close()
    output.close()


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
    # 2) Merge map data               #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Merge map data")
    merge_map_data(config["pcall"], config["map"], config["output"])

    print(">> Done.")

