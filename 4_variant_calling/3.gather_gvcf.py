# -*- coding: utf-8 -*-
"""
This script gathers the genotypes per chromosomes into one file
    then combines them back per scaffold/chromosomes into one file
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
from variable import *
from ntfy import Generate_cURL_cmd
import pandas as pd

# Directories:
direct = "{}/{}/vcf_files/".format(path, sp)
scaff_dir = "{}/{}/".format(path, sp)

# Import chrom/scaff names:
genome_assembly = pd.read_csv('{}scaffolds.txt'.format(scaff_dir),sep=' ', index_col=None, header=None)

# List of the files:
back_dir = []
geno_dir = []
for line in range(0,nb_scaff):
    scaff=genome_assembly.loc[line,1]
    back_dir.append("back_combine_genomicDBI_{}.g.vcf".format(scaff))
    geno_dir.append("genotype_genomicDBI_{}.g.vcf".format(scaff))


# The function:
def gather(direct, list_file, output, what, scaff):
    """Gather each chrom/scaff together"""
    gather_cmd = "gatk --java-options \"-XX:ParallelGCThreads=16 -Xmx120g \" GatherVcfs "
    for i in list_file:
        gather_cmd += "-I {}{} ".format(direct, i)
    gather_cmd += "-O {} ".format(output)
    """Create a .sh files with the gather functions."""
    file = open('{}gather_genomicDBImport_{}.sh'.format(direct, what),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 124G \n')
    file.write('#SBATCH --cpus-per-task=16 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write(Generate_cURL_cmd("Gathering {} {}".format(what, scaff)))
    file.write(gather_cmd)
    file.write(Generate_cURL_cmd("Finished gathering {} {}\nCommand: {}\nResult: $?".format(what, scaff, gather_cmd)))
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}gather_genomicDBImport_{}.out {}gather_genomicDBImport_{}.sh".format(direct, what, direct, what)
    subprocess.call(sub_cmd, shell=True)


################
# What you run #
################

# Gather if everything exist for all scaffolds/chromosomes:
list_exist=[]
for line in range(0, nb_scaff):
    scaff=genome_assembly.loc[line,1]
    list_exist.append(os.path.exists("{}back_combine_genomicDBI_{}.g.vcf".format(direct, scaff)))
if all(list_exist):
    print("\t All back_combine directories exist --> back combine done")
    mv_bc= "mv {}back_combine_genomicDBImport_* {}back_com.log/".format(direct, direct)
    subprocess.call(mv_bc, shell=True)
    print("\t Move the back combine log files")
    gather(direct=direct, list_file=back_dir, output="{}back_combine_genomicDBI_gather.g.vcf".format(direct), what="back_combine", scaff=scaff)
    print("Gather back combine")

list_exist=[]
for line in range(0, nb_scaff):
    scaff=genome_assembly.loc[line,1]
    list_exist.append(os.path.exists("{}genotype_genomicDBI_{}.g.vcf".format(direct, scaff)))
if all(list_exist):
    print("\t All genotypes directories exist --> genotype done")
    mv_geno= "mv {}genotype_genomicDBImport_* {}back_com.log/".format(direct, direct)
    subprocess.call(mv_geno, shell=True)
    print("\t Move the genotype log files")
    gather(direct=direct, list_file=geno_dir, output="{}genotype_genomicDBI_gather.g.vcf".format(direct), what="genotype", scaff=scaff)
    print("Gather genotypes")
