
# -*- coding: utf-8 -*-
"""
This script can be run when back has been splitted per trios and:
       - move the previous log files
       - find the mean depth save it in a file
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
from variable import *

# Directories:
direct = "{}/{}/vcf_handling/".format(path, sp)
ref="{}/{}/ref_fasta/{}.fa".format(path, sp,refGenome)
dir_tab="{}/{}/".format(path, sp)
recomb = "{}/{}/recombination/".format(path, sp)

# Dictionary:
f = open('{}/{}/pedigree.ped'.format(path, sp))
trio_dir = {}
for line in f:
    off = line.split()[1]
    fa = line.split()[2]
    mo = line.split()[3]
    name = off
    if name not in trio_dir:
        trio_dir[name] = []
    trio_dir[name].append((off, fa, mo))

def filter_trio(ref, off, fa, mo, direct, output_trio):
    """Select a trio"""
    filt_cmd = "gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx15g \" SelectVariants "
    filt_cmd += "-R {} ".format(ref)
    filt_cmd += "-V {}genotype_genomicDBI_{}_snp_filt.g.vcf ".format(direct, name)
    filt_cmd += "--exclude-filtered TRUE "
    filt_cmd += "--exclude-intervals {}genotype_genomicDBI_{}_snp_filt_MV.g.vcf ".format(direct, name)
    filt_cmd += "-O {}trio_{}_final_filtered_noMV.vcf.gz ".format(recomb, name)
    """Create a .sh files with the filter trio function."""
    file = open('{}{}.sh'.format(recomb, output_trio),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 15G \n')
    file.write('#SBATCH --cpus-per-task=1 \n')
    file.write('#SBATCH --time=2:59:00 \n')
    file.write(filt_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}{}.out {}{}.sh".format(recomb, output_trio, recomb, output_trio)
    subprocess.call(sub_cmd, shell=True)

for name in trio_dir:
    off = trio_dir[name][0][0]
    fa = trio_dir[name][0][1]
    mo = trio_dir[name][0][2]
    filter_trio(ref=ref, off= off, fa=fa, mo=mo, direct=direct, output_trio="filter_trio_{}".format(name))

print("It's been done")
