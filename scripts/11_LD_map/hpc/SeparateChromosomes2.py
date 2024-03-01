#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# SeparateChromosomes2.py
# -----------------------
# Run Lep-MAP3 SeparateChromosomes2 on Puhti.
# (HPC SCRIPT --> run wrapper)
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
    parser.add_argument("--scratch", "-scratch", help="Scratch path")
    parser.add_argument("--lepmap3", "-lepmap3", help="Lep-MAP3 path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--lodlimit", "-lodlimit", help="LOD limit")
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

### Import the datasets ###
def import_datasets( scratch_path, population, version ):
    param_file = "parameters_"+population+"_"+version+".txt"
    pcall_file = population+"_"+version+"_parent_call_filtered.txt"
    os.system("cp "+scratch_path+"/"+param_file+" "+param_file)
    os.system("cp "+scratch_path+"/"+pcall_file+" "+pcall_file)
    assert_creation(param_file)
    assert_creation(pcall_file)

### Load the list of parameter values ###
def load_parameters( population, version ):
    filename   = "parameters_"+population+"_"+version+".txt"
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
def run_LepMAP3_SeparateChromosomes2( scratch_path, lepmap3_path, population, version, lodlimit, parameters ):
    data_file   = population+"_"+version+"_parent_call_filtered.txt"
    output_file = population+"_"+version+"_LOD"+str(lodlimit)+"_map.txt"
    assert_creation(data_file)
    cmdline = []
    cmdline.append("java -cp "+lepmap3_path+" SeparateChromosomes2")
    cmdline.append("data="+data_file)
    cmdline.append("lodLimit="+str(lodlimit))
    cmdline.append("distortionLod="+parameters["distortionLod"])
    cmdline.append("numThreads="+parameters["numThreads"])
    cmdline.append("> "+output_file)
    os.system(" ".join(cmdline))
    assert_creation(output_file)
    os.system("gzip -c "+output_file+" > "+output_file+".gz")
    assert_creation(output_file+".gz")
    os.system("mv "+output_file+".gz "+scratch_path+"/"+output_file+".gz")
    assert_creation(scratch_path+"/"+output_file+".gz")


##################
#      MAIN      #
##################

if __name__ == '__main__':
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import the dataset           #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import datasets")
    import_datasets(config["scratch"], config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Load mapping parameters      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load mapping parameters")
    parameters = load_parameters(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Separate chromosomes         #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Separate chromosomes")
    run_LepMAP3_SeparateChromosomes2(config["scratch"], config["lepmap3"], config["population"], config["version"], config["lodlimit"], parameters)

    print(">> Done.")

