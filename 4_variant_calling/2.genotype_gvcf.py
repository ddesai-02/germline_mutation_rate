# -*- coding: utf-8 -*-
"""
This script genotypes and constructs back combined files
from the GenomicDBI folders VCF files per chromosome or scaffolds for all individuals.
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
ref = "{}/{}/ref_fasta/{}.fa".format(path, sp, refGenome)
scaff_dir = "{}/{}/".format(path, sp)

# Import chrom/scaff names:
genome_assembly = pd.read_csv('{}scaffolds.txt'.format(scaff_dir),sep=' ', index_col=None, header=None)

# The function:
def genotype(ref, sp, scaff, direct):
    """Genotypes all samples from GenomicsDBImport"""
    gt_cmd = "gatk --java-options \"-XX:ParallelGCThreads=3 -Xmx60g \" GenotypeGVCFs "
    gt_cmd += "-R {} ".format(ref)
    gt_cmd += "-V gendb://{}genomicDBI_{} ".format(scratch_dir, scaff)
    gt_cmd += "-G StandardAnnotation -new-qual "
    gt_cmd += "-O {}genotype_genomicDBI_{}.g.vcf ".format(direct, scaff)
    """Create a .sh files with the genotype functions."""
    file = open('{}genotype_genomicDBImport_{}.sh'.format(direct, scaff),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 64G \n')
    file.write('#SBATCH --cpus-per-task=3 \n')
    file.write('#SBATCH --time=12:00:00 \n')
    file.write('cp -a {}genomicDBI_{} {}genomicDBI_{} \n'.format(direct, scaff, scratch_dir, scaff))
    file.write(gt_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}genotype_genomicDBImport_{}.out {}genotype_genomicDBImport_{}.sh".format(direct, scaff, direct, scaff)
    subprocess.call(sub_cmd, shell=True)


def back_combine(ref, sp, scaff, direct):
    """Combine readable from GenomicsDBImport"""
    bc_cmd = "gatk --java-options \"-XX:ParallelGCThreads=3 -Xmx150g \" SelectVariants "
    bc_cmd += "-R {} ".format(ref)
    bc_cmd += "-V gendb://{}genomicDBI_{} ".format(scratch_dir, scaff)
    bc_cmd += "-O {}back_combine_genomicDBI_{}.g.vcf ".format(direct, scaff)
    """Create a .sh files with the genotype functions."""
    file = open('{}back_combine_genomicDBImport_{}.sh'.format(direct, scaff),'w')
    file.write('#!/bin/bash \n')
    file.write('#SBATCH --account={} \n'.format(account))
    file.write('#SBATCH --mem 156G \n')
    file.write('#SBATCH --cpus-per-task=3 \n')
    file.write('#SBATCH --time=10:00:00 \n')
    file.write('cp -a {}genomicDBI_{} /scratch/$SLURM_JOBID/genomicDBI_{} \n'.format(direct, scaff, scaff))
    file.write(bc_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .sh to the server"""
    sub_cmd = "sbatch -o {}back_combine_genomicDBImport_{}.out {}back_combine_genomicDBImport_{}.sh".format(direct, scaff, direct, scaff)
    subprocess.call(sub_cmd, shell=True)

################
# What you run #
################

# For each chromosome/scaffold one function:
list_exist=[]
for line in range(0, nb_scaff):
    scaff=genome_assembly.loc[line,1]
    list_exist.append(os.path.exists("{}genomicDBI_{}".format(direct, scaff)))
if all(list_exist):
    print("\t All the genomicDBI directoris exist --> combine variant done")
    mv_com = "mv {}combine_genomicDBImport_* {}combine.log".format(direct, direct)
    subprocess.call(mv_com, shell=True)
    print("\t Move the combine log files")
    mv_vcf = "mv {}*_res.g.vcf* {}inter_vcf".format(direct, direct)
    subprocess.call(mv_vcf, shell=True)
    print("\t Move the inter vcf files")
    for line in range(0, nb_scaff):
        scaff=genome_assembly.loc[line,1]
        genotype(ref=ref, sp=sp, scaff=scaff, direct=direct)
        print("Genotyping for {}".format(scaff))
        back_combine(ref=ref, sp=sp, scaff=scaff, direct=direct)
        print("Back combining to have combine for {} readable".format(scaff))

