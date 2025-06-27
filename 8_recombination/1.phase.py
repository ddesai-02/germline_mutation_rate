
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
recomb="{}/{}/recombination/".format(path, sp)
bam_dir="{}/{}/bam_files/".format(path, sp)

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

def phase_trio(ref, off, fa, mo, direct, output_trio):
    """Select a trio"""
    phase_cmd = "whatshap phase "
    phase_cmd += "--ped {}trio_{}.ped ".format(recomb, name)
    phase_cmd += "--reference {} ".format(ref)
    phase_cmd += "--output {}trio_{}_phased.vcf.gz ".format(recomb, name)
    phase_cmd += "{}trio_{}_final_filtered_noMV.vcf.gz ".format(recomb, name)
    phase_cmd += "{}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam {}{}_sorted.merged.addg.uniq.rmdup.bam".format(bam_dir, off, bam_dir, fa, bam_dir, mo)
    """Create a .sh files with the filter trio function."""
    file = open('{}{}.sh'.format(recomb, output_trio),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 15G \n')
    file.write('#SBATCH --cpus-per-task=8 \n')
    file.write('#SBATCH --time=11:59:00 \n')
    file.write(phase_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}{}.out {}{}.sh".format(recomb, output_trio, recomb, output_trio)
    subprocess.call(sub_cmd, shell=True)
    mv_cmd = "mv {}filter* {}filt.log".format(recomb,recomb)
    subprocess.call(mv_cmd, shell=True)

for name in trio_dir:
    off = trio_dir[name][0][0]
    fa = trio_dir[name][0][1]
    mo = trio_dir[name][0][2]
    file = open('{}trio_{}.ped'.format(recomb, name),'w')
    file.write('FAMILY01 {} {} {} -9 -9'.format(off, fa, mo))
    file.close()

for name in trio_dir:
    off = trio_dir[name][0][0]
    fa = trio_dir[name][0][1]
    mo = trio_dir[name][0][2]
    phase_trio(ref=ref, off= off, fa=fa, mo=mo, direct=direct, output_trio="phase_{}".format(name))

print("It's been done")
