#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2024 Charles Rocabert, Frédéric Guillaume
# Github: charlesrocabert/Tribolium-castaneum-transcriptomics-pipeline
#
# 2_RunParentCallFiltering.py
# ---------------------------
# Run Lep-MAP3 ParentCall2 and Filtering2 functions.
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
    parser.add_argument("--lepmap3", "-lepmap3", help="LepMap3 path")
    parser.add_argument("--population", "-population", help="Sample list population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--suffix", "-suffix", help="Suffix")
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

### Run Lep-MAP3 ParentCall2 function ###
def run_LepMAP3_ParentCall2( lepmap3_path, population, version, suffix, parameters ):
    pedigree_file = "./data/tribolium_pedigree/pedigree_"+population+"_"+version+".txt"
    vcf_file      = "./data/tribolium_snp/Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf"
    output_file   = "./data/tribolium_ld/"+population+"_"+version+"_parent_call.txt"
    cmdline       = []
    cmdline.append("java -cp "+lepmap3_path+" ParentCall2")
    cmdline.append("data="+pedigree_file)
    cmdline.append("vcfFile="+vcf_file)
    cmdline.append("halfSibs="+parameters["ParentCall2_halfSibs"])
    cmdline.append("removeNonInformative="+parameters["removeNonInformative"])
    cmdline.append("> "+output_file)
    os.system(" ".join(cmdline))
    assert_creation(output_file)
    os.system("gzip -c "+output_file+" > "+output_file+".gz")

### Run Lep-MAP3 Filtering2 function ###
### dataTolerance=0.01 suitable for multi-family data (official doc)
def run_LepMAP3_Filtering2( lepmap3_path, population, version, parameters ):
    data_file   = "./data/tribolium_ld/"+population+"_"+version+"_parent_call.txt"
    output_file = "./data/tribolium_ld/"+population+"_"+version+"_parent_call_filtered.txt"
    assert_creation(data_file)
    cmdline = []
    cmdline.append("java -cp "+lepmap3_path+" Filtering2")
    cmdline.append("data="+data_file)
    cmdline.append("removeNonInformative="+parameters["removeNonInformative"])
    cmdline.append("dataTolerance="+parameters["Filtering2_dataTolerance"])
    cmdline.append("> "+output_file)
    os.system(" ".join(cmdline))
    assert_creation(output_file)
    os.system("gzip -c "+output_file+" > "+output_file+".gz")


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
    # 2) Load mapping parameters      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Load mapping parameters")
    parameters = load_parameters(config["population"], config["version"])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Parent call                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parent call")
    run_LepMAP3_ParentCall2(config["lepmap3"], config["population"], config["version"], config["suffix"], parameters)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 4) Filtering                    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Filtering")
    run_LepMAP3_Filtering2(config["lepmap3"], config["population"], config["version"], parameters)

    print(">> Done.")

