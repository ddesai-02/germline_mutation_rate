"""
This script calculates the bases covered with depth >/= 5 of each individual bam file
and retrieves the length of the reference genome.
"""

import subprocess
import os
from variable import *

# Directories:
directory = "{}/{}/bam_files/".format(path, sp)
ref_dir = "{}/{}/ref_fasta/".format(path, sp)

# Retrieve names
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
name_dir = {}
for line in f:
		name = line.split()[0]
		if name not in name_dir:
			name_dir[name] = []
		name_dir[name] = "{}/{}/bam_files/{}_sorted.merged.addg.uniq.rmdup.bam".format(path, sp, name)


def breadth(input_bam, direct):
	cov_cmd = "samtools mpileup {} | awk -v X=\"$5\" '$4>=X' | wc -l".format(input_bam)
	file = open('{}_bases_cov.sh'.format(direct),'w')
	file.write('#!/bin/bash \n')
	file.write('#SBATCH --account={} \n'.format(account))
	file.write('#SBATCH --mem 8G \n')
	file.write('#SBATCH --cpus-per-task=2 \n')
	file.write('#SBATCH --time=00:50:00 \n')
	file.write(cov_cmd)
	file.write('\n')
	file.close()
	#""Submit the .sh to the server ""
	sub_cmd = "sbatch -o {}_breadth.out {}_bases_cov.sh".format(direct, direct)
	subprocess.call(sub_cmd, shell=True)

def ref_length(ref_dir, refGenome, directory):
	# Index Reference
	idx_cmd = "bowtie2-build {}{}.fa {}{}refgenome".format(ref_dir, refGenome, ref_dir, refGenome)
	# Reference Length
	len_cmd = "bowtie2-inspect -s {}refgenome | awk '{{ FS = \"\\t\" }} ; BEGIN{{L=0}}; {{L=L+$3}}; END{{print L}}'".format(ref_dir)
	file = open('{}ref_length.sh'.format(directory),'w')
	file.write('#!/bin/bash \n')
	file.write('#SBATCH --account={} \n'.format(account))
	file.write('#SBATCH --mem 8G \n')
	file.write('#SBATCH --cpus-per-task=8 \n')
	file.write('#SBATCH --time=02:59:00 \n')
	file.write(idx_cmd)
	file.write('\n')
	file.write(len_cmd)
	file.write('\n')
	file.close()
    #""Submit the .sh to the server ""
	sub_cmd = "sbatch -o {}ref_length.out {}ref_length.sh".format(directory, directory)
	subprocess.call(sub_cmd, shell=True)

for name in name_dir:
	breadth(input_bam="{}{}_sorted.merged.addg.uniq.rmdup.bam".format(directory, name),direct="{}{}".format(directory, name))

ref_length(ref_dir=ref_dir, refGenome=refGenome, directory=directory)
