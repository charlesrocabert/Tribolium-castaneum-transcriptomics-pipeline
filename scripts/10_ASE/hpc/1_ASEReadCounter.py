#!/usr/bin/env python3
# coding: utf-8

#***************************************************************************
# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume
# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation
#
# 1_ASEReadCounter.py
# -------------------
# Run GATK ASEReadCounter for a given sample.
# (HPC SCRIPT --> array wrapper)
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
    parser.add_argument("--bucket", "-bucket", help="Allas bucket name")
    parser.add_argument("--population", "-population", help="Population")
    parser.add_argument("--version", "-version", help="Reference genome version")
    parser.add_argument("--suffix", "-suffix", help="Suffix")
    parser.add_argument("--sindex", "-sindex", type=int, default=0, help="Index of the sample in the list of samples")
    args = parser.parse_args()
    return(vars(args))

### Check the creation of a file ###
def assert_creation( filename ):
    assert os.path.isfile(filename), ">> "+filename+" has not been created. Exit."

### Check the deletion of a file ###
def assert_deletion( filename ):
    assert not os.path.isfile(filename), ">> "+filename+" has not been deleted. Exit."

### Import VCF file ###
def import_sample_list( population, version ):
    filename = "samples_"+population+"_"+version+".csv"
    os.system("cp /scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_bam/"+filename+" .")
    assert_creation(filename)

### Load the list of samples ###
def load_sample_list( population, version ):
    samples   = []
    file      = open("samples_"+population+"_"+version+".csv", "r")
    csvreader = csv.reader(file, delimiter=";")
    header    = next(csvreader)
    samples.append(header)
    for row in csvreader:
        sample = {}
        for i in range(len(header)):
            sample[header[i]] = row[i]
        samples.append(sample)
    file.close()
    return samples

### Import BAM file ###
def import_bam( bucket, sample_name ):
    filename = sample_name+".bam"
    os.system("a-get "+bucket+"/bam/"+filename)
    assert_creation(filename)

### Create BAM file index ###
def create_bam_index( sample_name ):
    os.system("samtools index "+sample_name+".bam")
    assert_creation(sample_name+".bam.bai")

### Import VCF file ###
def import_VCF( population, version, suffix ):
    filename = "Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf.gz"
    os.system("cp /scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_snp/"+filename+" .")
    assert_creation(filename)

### Create the VCF index file ###
def create_VCF_index( population, version, suffix ):
    vcf_file = "Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf.gz"
    os.system("gunzip "+vcf_file)
    os.system("bgzip -c "+vcf_file.strip(".gz")+" > "+vcf_file)
    os.system("tabix -p vcf "+vcf_file)
    os.system("gatk IndexFeatureFile -I "+vcf_file)
    assert_creation(vcf_file+".tbi")
    #assert_creation(vcf_file+".idx")

### Import reference genome ###
def import_reference_genome( version ):
    filename = "Tribolium_castaneum_"+version+".fna"
    os.system("cp /scratch/project_2003847/Tribolium-Polygenic-Adaptation/data/tribolium_genome/Tribolium_castaneum_"+version+"/"+filename+" .")
    assert_creation(filename)

### Create reference genome indices ###
def create_reference_genome_indices( version ):
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    os.system("samtools faidx "+reference_genome)
    os.system("gatk CreateSequenceDictionary -R "+reference_genome)
    assert_creation(reference_genome+".fai")
    assert_creation(reference_genome.strip(".fna")+".dict")

### Run GATK ASEReadCounter ###
def run_GATK_ASEReadCounter( population, version, suffix, sample_name ):
    reference_genome = "Tribolium_castaneum_"+version+".fna"
    vcf_file         = "Tribolium_castaneum_"+population+"_"+version+"_"+suffix+".vcf.gz"
    bam_file         = sample_name+".bam"
    output_file      = population+"_"+version+"_"+suffix+"_"+sample_name+".table"
    cmdline          = []
    cmdline.append("gatk ASEReadCounter")
    cmdline.append("-R "+reference_genome)
    cmdline.append("-V "+vcf_file)
    cmdline.append("-I "+bam_file)
    cmdline.append("-O "+output_file)
    os.system(" ".join(cmdline))
    assert_creation(output_file)

### Export count table to Allas ###
def export_table( population, version, suffix, sample_name ):
    output_file = population+"_"+version+"_"+suffix+"_"+sample_name+".table"
    os.system("cp "+output_file+" /scratch/project_2003847/Tribolium_castaneum_ASE/tables/"+output_file)


##################
#      MAIN      #
##################

if __name__ == '__main__':
    print("")
    print("#***************************************************************************")
    print("# Copyright © 2021-2023 Charles Rocabert, Frédéric Guillaume")
    print("# Web: https://github.com/charlesrocabert/Tribolium-Polygenic-Adaptation")
    print("#")
    print("# 1_ASEReadCounter.py")
    print("# -------------------")
    print("# Run GATK ASEReadCounter for a given sample.")
    print("# (HPC SCRIPT --> array wrapper)")
    print("#***************************************************************************")
    print("")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1) Parse command line arguments                  #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Parse command line arguments")
    config = parse_arguments()
    os.chdir(os.environ['DATADIR'])

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 2) Import the list of samples and get the sample #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    print(">> Import the list of samples")
    import_sample_list(config["population"], config["version"])
    samples     = load_sample_list(config["population"], config["version"])
    sample      = samples[config["sindex"]]
    sample_name = sample["sample"]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 3) Run GATK ASEReadCounter                       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    ### 3.1) Import the BAM file ###
    print(">> Import BAM file and create indices")
    import_bam(config["bucket"], sample_name)
    create_bam_index(sample_name)

    ### 3.2) Import the VCF file ###
    print(">> Import VCF index and create indices")
    import_VCF(config["population"], config["version"], config["suffix"])
    create_VCF_index(config["population"], config["version"], config["suffix"])

    ### 3.3) Import the reference genome and create index ###
    print(">> Import reference genome and create indices")
    import_reference_genome(config["version"])
    create_reference_genome_indices(config["version"])

    ### 3.4) Run GATK ASEReadCounter ###
    print(">> Run GATK ASEReadCounter")
    run_GATK_ASEReadCounter(config["population"], config["version"], config["suffix"], sample_name)

    ### 3.5) Export the resulting count table ###
    print(">> Export count table")
    export_table(config["population"], config["version"], config["suffix"], sample_name)

    print(">> Done.")

