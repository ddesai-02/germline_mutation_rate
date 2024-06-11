# -*- coding: utf-8 -*-
"""
This script call variants in BP RESOLUTION for all individuals per chromosome/scaffold
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
ref_dir = "{}/{}/ref_fasta/{}.fa".format(path, sp, refGenome)
bam_dir = "{}/{}/bam_files/".format(path, sp)
vcf_dir = "{}/{}/vcf_files/".format(path, sp)
scaff_dir = "{}/{}/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam".format(path, sp, name)


# Import chrom or scaffold names:
genome_assembly = pd.read_csv('{}scaffolds.txt'.format(scaff_dir),sep=' ', index_col=None, header=None)


# The function:
def call_var(ref, in_bam, out_vcf, scaff_old, out_dir):
    """Haplotype caller function to call variants for each samples"""
    call_cmd = "gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx90g -Djava.io.tmpdir={}\" HaplotypeCaller ".format(scratch_dir)
    call_cmd += "-R {} ".format(ref)
    call_cmd += "-I {} ".format(in_bam)
    call_cmd += "-O {} ".format(out_vcf)
    call_cmd += "-ERC BP_RESOLUTION "
    call_cmd += "-L {} ".format(scaff_old)
    call_cmd += "--dont-use-soft-clipped-bases "
    call_cmd += "--native-pair-hmm-threads 1 "
    call_cmd += "--tmp-dir {} ".format(scratch_dir)
    """Create a .sh files with the calling variant functions."""
    file = open('{}_call_res_g_{}.sh'.format(out_dir, scaff_old),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account{}} \n'.format(account))
    file.write('#SBATCH --mem 16G \n')
    file.write('#SBATCH --cpus-per-task=1 \n')
    file.write('#SBATCH --time=23:00:00 \n')
    file.write(call_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}_call_res_g_{}.out {}_call_res_g_{}.sh".format(out_dir, scaff_old, out_dir, scaff_old)
    subprocess.call(sub_cmd, shell=True)


#################
# What you run  #
#################

vcf_files_dir = open("{}/{}/vcf_files.txt".format(path, sp), "w")

for name in bamfile_dir: # for each individual
    for scaff in range(0, len(genome_assembly[0]) - 1): # for each scaffold
        current = "{}_{}".format(name, genome_assembly[0][scaff])
        print(current)

        if os.path.exists("{}{}_res.g.vcf".format(vcf_dir, current)): # check if vcf for individual exists
            print("\t The res.g.vcf file for {} already exists --> CALL VARIANT DONE".format(current))
        else:
            print("\t The res.g.vcf file for {} doesn't exist --> submit the function".format(current))
            scaff_old=genome_assembly[0][scaff]
            call_var(ref=ref_dir, in_bam=bamfile_dir[name], out_vcf="{}{}_{}_res.g.vcf".format(vcf_dir, name, scaff_old), scaff_old=scaff_old, out_dir="{}{}".format(vcf_dir, name))
            vcf_files_dir.write(current + "_res.g.vcf \n")
vcf_files_dir.close()

