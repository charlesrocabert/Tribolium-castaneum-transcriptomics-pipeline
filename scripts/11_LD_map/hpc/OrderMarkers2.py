#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# OrderMarkers2.py
# ----------------
# Run Lep-MAP3 OrderMarkers2 on Puhti.
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
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--in-pcall", "-in-pcall", help="Input parent call suffix")
    parser.add_argument("--in-map", "-in-map", help="Input map suffix")
    parser.add_argument("--out-map", "-out-map", help="Output map suffix")
    parser.add_argument("--chr", "-chr", help="Chromosome to order")
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
def import_datasets( scratch_path, population, version, in_pcall, in_map ):
    param_file = "parameters_"+population+"_"+version+".txt"
    pcall_file = population+"_"+version+"_"+in_pcall+".txt"
    map_file   = population+"_"+version+"_"+in_map+".txt"
    os.system("cp "+scratch_path+"/"+param_file+" "+param_file)
    os.system("cp "+scratch_path+"/"+pcall_file+" "+pcall_file)
    os.system("cp "+scratch_path+"/"+map_file+" "+map_file)
    assert_creation(param_file)
    assert_creation(pcall_file)
    assert_creation(map_file)

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

### Run Lep-MAP3 OrderMarkers2 function ###
def run_LepMAP3_OrderMarkers2( scratch_path, lepmap3_path, population, version, in_pcall, in_map, out_map, chr, parameters ):
    data_file   = population+"_"+version+"_"+in_pcall+".txt"
    map_file    = population+"_"+version+"_"+in_map+".txt"
    output_file = population+"_"+version+"_"+out_map+".txt"
    assert_creation(data_file)
    assert_creation(map_file)
    cmdline = []
    cmdline.append("java -cp "+lepmap3_path+" OrderMarkers2")
    cmdline.append("data="+data_file)
    cmdline.append("map="+map_file)
    cmdline.append("chromosome="+chr)
    cmdline.append("outputPhasedData="+parameters["OrderMarkers2_outputPhasedData"])
    cmdline.append("sexAveraged="+parameters["OrderMarkers2_sexAveraged"])
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
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# OrderMarkers2.py")
    print("# ----------------")
    print("# Run Lep-MAP3 OrderMarkers2 on Puhti.")
    print("# (HPC SCRIPT --> run wrapper)")
    print("#***************************************************************************")
    print("")

    SCRATCH_PATH = "/scratch/project_2003847/Tribolium_castaneum_ldMAP"
    LEPMAP3_PATH = "/scratch/project_2003847/Tribolium_castaneum_ldMAP/Lep-MAP3/bin"
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
    import_datasets(SCRATCH_PATH, config["population"], config["version"], config["in_pcall"], config["in_map"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Load mapping parameters      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load mapping parameters")
    parameters = load_parameters(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Order markers                #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Order markers")
    run_LepMAP3_OrderMarkers2(SCRATCH_PATH, LEPMAP3_PATH, config["population"], config["version"], config["in_pcall"], config["in_map"], config["out_map"], config["chr"], parameters)

    print(">> Done.")

