#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 3_RunSeparateChromosomes.py
# ---------------------------
# Run Lep-MAP3 SeparateChromosomes2.
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
    parser.add_argument("--out-map", "-out-map", help="Output map suffix")
    parser.add_argument("--lodlimit", "-lodlimit", help="LOD limit value")
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

### Run Lep-MAP3 SeparateChromosomes2 function ###
def run_LepMAP3_SeparateChromosomes2( lepmap3_path, population, version, in_pcall, out_map, lodlimit, parameters ):
    data_file   = "./data/tribolium_ld/"+population+"_"+version+"_"+in_pcall+".txt"
    output_file = "./data/tribolium_ld/"+population+"_"+version+"_"+out_map+".txt"
    assert_creation(data_file)
    cmdline = []
    cmdline.append("java -cp "+lepmap3_path+" SeparateChromosomes2")
    cmdline.append("data="+data_file)
    #cmdline.append("lodLimit="+parameters["SeparateChromosomes2_lodLimit"])
    cmdline.append("lodLimit="+lodlimit)
    #cmdline.append("minLod="+parameters["SeparateChromosomes2_minLod"]) # Not useful
    cmdline.append("distortionLod="+parameters["distortionLod"])
    #cmdline.append("phasedData="+parameters["SeparateChromosomes2_phasedData"]) # Useful only if you have grand parents
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
    print("# 3_RunSeparateChromosomes.py")
    print("# ---------------------------")
    print("# Run Lep-MAP3 SeparateChromosomes2.")
    print("# (LOCAL SCRIPT)")
    print("#***************************************************************************")
    print("")

    WD_PATH      = "/Users/charlesrocabert/git/Tribolium-Polygenic-Adaptation"
    LEPMAP3_PATH = "/Users/charlesrocabert/Lep-MAP3/bin"
    os.chdir(WD_PATH)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Load mapping parameters      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load mapping parameters")
    parameters = load_parameters(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Separate chromosomes         #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Separate chromosomes")
    run_LepMAP3_SeparateChromosomes2(LEPMAP3_PATH, config["population"], config["version"], config["in_pcall"], config["out_map"], config["lodlimit"], parameters)

    print(">> Done.")

