#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 5_RunJoinSingles2All.py
# -----------------------
# Run Lep-MAP3 JoinSingles2All.
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

### Run Lep-MAP3 JoinSingles2All function ###
def run_LepMAP3_JoinSingles2All( lepmap3_path, population, version, in_map, out_map, parameters ):
    data_file   = "./data/tribolium_ld/"+population+"_"+version+"_parent_call_filtered.txt"
    map_file    = "./data/tribolium_ld/"+population+"_"+version+"_"+in_map+".txt"
    output_file = "./data/tribolium_ld/"+population+"_"+version+"_"+out_map+".txt"
    assert_creation(data_file)
    assert_creation(map_file)
    cmdline = []
    cmdline.append("java -cp "+lepmap3_path+" JoinSingles2All")
    cmdline.append("data="+data_file)
    cmdline.append("map="+map_file)
    cmdline.append("lodLimit="+parameters["JoinSingles2All_lodLimit"])
    cmdline.append("distortionLod="+parameters["distortionLod"])
    cmdline.append("iterate="+parameters["JoinSingles2All_iterate"])
    cmdline.append("numThreads="+parameters["numThreads"])
    cmdline.append("> "+output_file)
    os.system(" ".join(cmdline))
    assert_creation(output_file)
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
    print("# 5_RunJoinSingles2All.py")
    print("# -----------------------")
    print("# Run Lep-MAP3 JoinSingles2All.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH      = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    LEPMAP3_PATH = "/Users/charlesrocabert/Lep-MAP3/bin"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments     #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load mapping parameters          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load mapping parameters")
    parameters = load_parameters(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Join single markers to LD groups #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Join single markers to LD groups")
    run_LepMAP3_JoinSingles2All(LEPMAP3_PATH, config["population"], config["version"], config["in_map"], config["out_map"], parameters)

    print(">> Done.")

