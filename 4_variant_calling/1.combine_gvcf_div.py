# -*- coding: utf-8 -*-
"""
This script catenate the GVCF files per scaffold/chromosome for all individuals.
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Directories:
direct = "{}/{}/vcf_files/".format(path, sp)
scaff_dir = "{}/{}/".format(path, sp)


# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam".format(path, sp, name)

# And VCF
f = open('{}/{}/vcf_files.txt'.format(path, sp))
vcf_dir = {}
for line in f:
    name = line.split()[0]
    if name not in vcf_dir:
        vcf_dir[name] = []
    vcf_dir[name] = "{}/{}/vcf_files/{}".format(path, sp, name)


# Import scaff names:
genome_assembly = pd.read_csv('{}scaffolds.txt'.format(scaff_dir),sep=' ', index_col=None, header=None)

# The function:
def combine(scaffold, direct):
    """Combine all samples with GenomicsDBImport"""
    combine_cmd = "gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx100g -Djava.io.tmpdir={}\" GenomicsDBImport ".format(scratch_dir)
    for i in vcf_dir:
        combine_cmd += "--variant {}{} ".format(direct,i)
    combine_cmd += "--tmp-dir {} ".format(scratch_dir)
    combine_cmd += "--genomicsdb-workspace-path {}genomicDBI_{} ".format(scratch_dir, scaffold)
    combine_cmd += "-L {} ".format(scaffold)
    """Create a .sh files with the combine variant functions."""
    file = open('{}combine_genomicDBImport_{}.sh'.format(direct, scaffold),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 110G \n')
    file.write('#SBATCH --cpus-per-task=1 \n')
    file.write('#SBATCH --time=10:00:00 \n')
##    file.write('#SBATCH --time=250:00:00 \n')
    file.write(combine_cmd)
    file.write('\n')
    file.write('cp -a {}genomicDBI_{} {}genomicDBI_{}'.format(scratch_dir, scaffold, direct, scaffold))
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}combine_genomicDBImport_{}.out {}combine_genomicDBImport_{}.sh".format(direct, scaffold, direct, scaffold)
    subprocess.call(sub_cmd, shell=True)


################
# What you run #
################

# For each scaffosome one function:
list_exist=[]
for file in list(vcf_dir.values()):
    list_exist.append(os.path.exists(file))
if all(list_exist):
    print("\t All the res.g.vcf files exist --> call variant done")
    mv_call = "mv {}*_call_res_g* {}call.log".format(direct, direct)
    subprocess.call(mv_call, shell=True)
    print("\t Move the call log files")
    for scaff in range(0, nb_scaff):
        scaffold=genome_assembly[0][scaff]
        vcf_dir=[]
        for name in bamfile_dir:
            vcf_dir.append("{}_{}_res.g.vcf".format(name, scaffold))
        print("Combine for scaffold/chromosome {} called {}".format(scaffold, scaffold))
        combine(scaffold=scaffold, direct=direct)
else:
    print("\t Some res.g.vcf file are missing --> PROBLEM")