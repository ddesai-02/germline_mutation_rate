"""
This script calculates the breadth of each individual bam file.
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


file = open("{}breadths_out.txt".format(directory), 'w')
for name in name_dir:
	covered = subprocess.check_output(['tail', '-1', '{}{}_breadth.out'.format(directory, name)])
	covered = covered.decode("utf-8").strip('\n')
	print('{}'.format(covered))
	print('{}'.format(name))
	file.write('{} {} \n'.format(name, covered))
file.close()

f = open('/home/devan/projects/rrg-shaferab/devan/Odocoileus_virginianus/bam_files/breadths_out.txt')
breadths = open('/home/devan/projects/rrg-shaferab/devan/Odocoileus_virginianus/breadths.txt', 'w')

for line in f:
	content = line.split()
	name = content[0]
	covered = int(content[1])
	ref_length = int(subprocess.check_output(['tail', '-1', '{}ref_length.out'.format(directory)]))
	breadth = covered/ref_length
	print(name, covered, breadth)
	breadths.write('{} {} {}\n'.format(name, covered, breadth))
