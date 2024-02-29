#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 6_CleanLinkageGroups.py
# ------------------------
# Filter out non informative linkage groups and SNPs.
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
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--in-pcall", "-in-pcall", help="Input parent call suffix")
    parser.add_argument("--in-map", "-in-map", help="Input map suffix")
    parser.add_argument("--out-pcall", "-out-pcall", help="Output parent call suffix")
    parser.add_argument("--out-map", "-out-map", help="Output map suffix")
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

### Load the list of parameter values ###
def load_parameters( population, version ):
    filename   = "./data/tribolium_ld/parameters_"+population+"_"+version+".txt"
    parameters = {}
    f          = open(filename, "r")
    l          = f.readline()
    while l:
        if not l.startswith("#") and len(l) >= 2:
            l = l.strip("\n").split(" ")
            assert l[0] not in parameters.keys()
            parameters[l[0]] = l[1]
        l = f.readline()
    f.close()
    return parameters

### Filter out non informative linkage groups ###
def clean_linkage_groups( population, version, in_pcall, in_map, out_pcall, out_map, LGs_to_collect ):
    parentcall_file         = "./data/tribolium_ld/"+population+"_"+version+"_"+in_pcall+".txt"
    map_file                = "./data/tribolium_ld/"+population+"_"+version+"_"+in_map+".txt"
    parentcall_cleaned_file = "./data/tribolium_ld/"+population+"_"+version+"_"+out_pcall+".txt"
    map_cleaned_file        = "./data/tribolium_ld/"+population+"_"+version+"_"+out_map+".txt"
    #----------------------------------------------------------#
    # 1) Collect the 10 first LGs with their chromosome counts #
    #----------------------------------------------------------#
    print("   > Collect linkage groups and chromosomes information")
    ### Open files ###
    LGs   = {}
    map   = open(map_file, "r")
    pcall = open(parentcall_file, "r")
    ### Manage headers ###
    for i in range(1):
        map.readline()
    for i in range(7):
        pcall.readline()
    ### Parse files and collect information ###
    l1    = map.readline()
    l2    = pcall.readline()
    count = 0
    while l1:
        assert l2
        if count % 10000 == 0:
            print("     • "+str(count)+" markers parsed")
        count += 1
        l1  = l1.strip("\n").split("\t")
        l2  = l2.strip("\n").split("\t")
        lg  = int(l1[0])
        chr = l2[0]
        if lg in LGs_to_collect:
            if lg not in LGs.keys():
                LGs[lg] = {}
            else:
                if chr not in LGs[lg].keys():
                    LGs[lg][chr] = 1
                else:
                    LGs[lg][chr] += 1
        l1 = map.readline()
        l2 = pcall.readline()
    assert not l2
    map.close()
    pcall.close()
    #----------------------------------------------------------#
    # 2) Identify the major chromosome per LG                  #
    #----------------------------------------------------------#
    print("   > Find major chromosome per LG")
    for lg in LGs_to_collect:
        assert lg in LGs.keys()
        best_count = 0
        best_chr   = ""
        for chr in LGs[lg].keys():
            if LGs[lg][chr] > best_count:
                best_count = LGs[lg][chr]
                best_chr   = chr
        LGs[lg]["major_chr"] = best_chr
    #----------------------------------------------------------#
    # 3) Rebuild the map file and the parent call file         #
    #----------------------------------------------------------#
    print("   > Rebuild files")
    pcall         = open(parentcall_file, "r")
    pcall_cleaned = open(parentcall_cleaned_file, "w")
    map           = open(map_file, "r")
    map_cleaned   = open(map_cleaned_file, "w")
    ### Manage headers ###
    for i in range(1):
        map_cleaned.write(map.readline())
    for i in range(7):
        pcall_cleaned.write(pcall.readline())
    ### Parse files ###
    l1    = map.readline()
    l2    = pcall.readline()
    count = 0
    while l1:
        assert l2
        if count % 10000 == 0:
            print("     • "+str(count)+" markers parsed")
        count += 1
        lg  = int(l1.strip("\n").split("\t")[0])
        chr = l2.strip("\n").split("\t")[0]
        if lg in LGs.keys():
            if chr == LGs[lg]["major_chr"]:
                map_cleaned.write(l1)
                pcall_cleaned.write(l2)
        l1 = map.readline()
        l2 = pcall.readline()
    assert not l2
    map.close()
    pcall.close()
    map_cleaned.close()
    pcall_cleaned.close()
    os.system("gzip -c "+parentcall_cleaned_file+" > "+parentcall_cleaned_file+".gz")
    os.system("gzip -c "+map_cleaned_file+" > "+map_cleaned_file+".gz")

##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 6_CleanLinkageGroups.py")
    print("# -----------------------")
    print("# Filter out non informative linkage groups.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH      = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    LEPMAP3_PATH = "/Users/charlesrocabert/Lep-MAP3/bin"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Clean linkage groups map file #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Clean linkage groups")
    clean_linkage_groups(config["population"], config["version"], config["in_pcall"], config["in_map"], config["out_pcall"], config["out_map"], range(1,11))

    print(">> Done.")

