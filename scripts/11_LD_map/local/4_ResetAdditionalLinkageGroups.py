#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 4_ResetAdditionalLinkageGroups.py
# ---------------------------------
# Filter out non informative linkage groups.
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
    parser.add_argument("--in-map", "-in-map", help="Input map suffix")
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

### Reset additional LGs to zero ###
def reset_additional_linkage_groups( population, version, in_map, out_map, LGs_to_keep ):
    map_file    = "./data/tribolium_ld/"+population+"_"+version+"_"+in_map+".txt"
    output_file = "./data/tribolium_ld/"+population+"_"+version+"_"+out_map+".txt"
    print("   > Reset additional linkage groups to zero")
    ### Open files ###
    map     = open(map_file, "r")
    new_map = open(output_file, "w")
    ### Manage headers ###
    for i in range(1):
        new_map.write(map.readline())
    ### Parse files and modify information ###
    l     = map.readline()
    count = 0
    while l:
        if count % 10000 == 0:
            print("     • "+str(count)+" markers parsed")
        count += 1
        l  = l.strip("\n").split("\t")
        LG = int(l[0])
        if LG in LGs_to_keep:
            new_map.write("\t".join(l)+"\n")
        else:
            l[0] = "0"
            new_map.write("\t".join(l)+"\n")
        l = map.readline()
    map.close()
    new_map.close()
    os.system("gzip -c "+output_file+" > "+output_file+".gz")

##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 4_ResetAdditionalLinkageGroups.py")
    print("# ---------------------------------")
    print("# Reset additional linkage groups to zero.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Clean linkage groups map file #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Reset additional linkage groups")
    reset_additional_linkage_groups(config["population"], config["version"], config["in_map"], config["out_map"], range(1,11))

    print(">> Done.")

