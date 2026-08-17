# -*- coding: utf-8 -*-
"""
This script creates all the files necessary to run the phasing for each trio.
It generates:
1. sites_{offspring_id}.txt: List of DNM sites (CHROM, POS) for the offspring.
2. dnm_windows_{offspring_id}.bed: BED file of genomic windows around each DNM.
3. dnm_windows_clipped_{offspring_id}.bed: Clipped BED file (windows adjusted to chromosome ends).
4. trio_genotypes_{offspring_id}.vcf.gz: VCF containing only the trio samples, bgzipped and indexed.
5. dnms_flank_{offspring_id}.vcf.gz: VCF file containing only DNMs within the defined windows for the trio.
Each set of operations for a trio is submitted as a separate SBATCH job.
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
import pandas as pd
from variable import * # Import your variables.py (sp, path, refGenome, account, scratch_dir)

# --- Directories ---
direct_denovo = "{}/{}/de_novo_mutation/".format(path, sp)
direct_phasing = "{}/{}/phasing/".format(path, sp)
ref_fasta_dir = "{}/{}/ref_fasta/".format(path,sp)
direct_vcf_files = "{}/{}/vcf_files/".format(path, sp) # Assuming main VCF is here
direct = "{}/{}/".format(path, sp) # Base project directory

# --- Parameters ---
PHASING_WINDOW_FLANK = 150 # Window size for flanking regions around DNMs (e.g., 150bp on each side)

# --- Dictionary for trios ---
pedigree_path = '{}pedigree.ped'.format(direct)
trio_data = {}
try:
    with open(pedigree_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4: # Skip malformed lines
                continue
            offspring_id = parts[1]
            mother_id = parts[3]
            father_id = parts[2]
            trio_data[offspring_id] = (offspring_id, mother_id, father_id)
except FileNotFoundError:
    print(f"Error: Pedigree file not found at {pedigree_path}. Please check the path in variable.py.")
    exit(1)


# --- The function to create and submit the prerequisite job for a single trio ---
def prepare_trio_phasing_files_sbatch(offspring_id, mother_id, father_id, 
                                      denovo_dir, phasing_dir, ref_fasta_dir, vcf_files_dir,
                                      window_flank, account_name, scratch_dir, ref_genome_name):
    """
    Creates a .sh script to generate prerequisite files for phasing a single trio and submits it via sbatch.
    """
    
    print(f"--- Preparing SBATCH job for prerequisite files for trio: {offspring_id} ---")

    # Ensure the output directory for this trio exists within the phasing directory
    trio_output_dir = os.path.join(phasing_dir, offspring_id)
    os.makedirs(trio_output_dir, exist_ok=True)

    # --- Input file paths ---
    dnm_tab_file = os.path.join(denovo_dir, f"data_denovo_{offspring_id}_new.tab")
    fasta_fai_file = os.path.join(ref_fasta_dir, f"{ref_genome_name}.fa.fai")
    main_vcf_file = os.path.join(vcf_files_dir, "genotype_genomicDBI_gather.g.vcf")

    # --- Output file paths ---
    sites_file = os.path.join(trio_output_dir, f"sites_{offspring_id}.txt")
    dnm_windows_bed = os.path.join(trio_output_dir, f"dnm_windows_{offspring_id}.bed")
    dnm_windows_clipped_bed = os.path.join(trio_output_dir, f"dnm_windows_clipped_{offspring_id}.bed")
    trio_vcf_gz = os.path.join(trio_output_dir, f"trio_genotypes_{offspring_id}.vcf.gz")
    dnms_flank_vcf_gz = os.path.join(trio_output_dir, f"dnms_flank_{offspring_id}.vcf.gz")
    
    sh_script_name = os.path.join(trio_output_dir, f"prepare_phasing_files_{offspring_id}.sh")
    log_file = os.path.join(trio_output_dir, f"prepare_phasing_files_{offspring_id}.log")
    err_file = os.path.join(trio_output_dir, f"prepare_phasing_files_{offspring_id}.err")

    # --- Create the .sh script ---
    with open(sh_script_name, 'w') as f_sh:
        f_sh.write('#!/bin/bash\n')
        f_sh.write(f'#SBATCH --account={account_name}\n')
        f_sh.write('#SBATCH --mem=8G\n') # Increased memory for VCF operations
        f_sh.write('#SBATCH --cpus-per-task=1\n') # Most of these commands are single-threaded
        f_sh.write('#SBATCH --time=01:00:00\n') # Increased time for potentially large VCFs
        f_sh.write(f'#SBATCH -o {log_file}\n')
        f_sh.write(f'#SBATCH -e {err_file}\n')
        f_sh.write('\n')
        f_sh.write(f'echo "Starting prerequisite file generation for trio: {offspring_id}"\n')
        f_sh.write(f'echo "Job ID: $SLURM_JOB_ID"\n')
        f_sh.write(f'echo "Date: $(date)"\n')
        f_sh.write('\n')
        
        # Add module loads if necessary for bcftools, bgzip, tabix on your cluster
        # e.g., f_sh.write('module load StdEnv/2023 bcftools\n')
        # Check with your cluster documentation or 'module avail'

        # 1. Create sites_{offspring_id}.txt
        f_sh.write(f'echo "Generating {os.path.basename(sites_file)}..."\n')
        f_sh.write(f"awk 'NR>1 {{print $2 \"\\t\" $3}}' {dnm_tab_file} > {sites_file}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: awk for sites_file failed."; exit 1; fi\n')
        f_sh.write('\n')

        # 2. Create dnm_windows_{offspring_id}.bed
        f_sh.write(f'echo "Generating {os.path.basename(dnm_windows_bed)}..."\n')
        f_sh.write(f"awk '{{start=($2-{window_flank}); if (start<0) start=0; print $1\"\\t\"start\"\\t\"($2+{window_flank})}}' {sites_file} > {dnm_windows_bed}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: awk for dnm_windows_bed failed."; exit 1; fi\n')
        f_sh.write('\n')

        # 3. Create dnm_windows_clipped_{offspring_id}.bed
        f_sh.write(f'echo "Generating {os.path.basename(dnm_windows_clipped_bed)}..."\n')
        f_sh.write(f"awk 'NR==FNR {{ maxlen[$1]=$2; next }} \\\n")
        f_sh.write(f"     {{ chrom=$1; start=$2; end=$3; \\\n")
        f_sh.write(f"       if (end > maxlen[chrom]) end = maxlen[chrom]; \\\n")
        f_sh.write(f"       print chrom \"\\t\" start \"\\t\" end; \\\n")
        f_sh.write(f"     }}' {fasta_fai_file} {dnm_windows_bed} > {dnm_windows_clipped_bed}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: awk for dnm_windows_clipped_bed failed."; exit 1; fi\n')
        f_sh.write('\n')
        
        # 4. Create trio-specific VCF and bgzip/index it
        f_sh.write(f'echo "Extracting, compressing, and indexing trio VCF: {os.path.basename(trio_vcf_gz)}..."\n')
        f_sh.write(f"bcftools view -s {mother_id},{father_id},{offspring_id} {main_vcf_file} -Oz -o {trio_vcf_gz}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: bcftools view for trio_vcf_gz failed."; exit 1; fi\n')
        f_sh.write(f"tabix -p vcf {trio_vcf_gz}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: tabix for trio_vcf_gz failed."; exit 1; fi\n')
        f_sh.write('\n')

        # 5. Create dnms_flank_{offspring_id}.vcf.gz (now using the trio-specific VCF)
        f_sh.write(f'echo "Generating {os.path.basename(dnms_flank_vcf_gz)}..."\n')
        f_sh.write(f"bcftools view -R {dnm_windows_clipped_bed} -Oz -o {dnms_flank_vcf_gz} {trio_vcf_gz}\n")
        f_sh.write('if [ $? -ne 0 ]; then echo "Error: bcftools view for dnms_flank_vcf_gz failed."; exit 1; fi\n')
        f_sh.write(f"tabix -p vcf {dnms_flank_vcf_gz}\n")
        f_sh.write('\n')
        f_sh.write(f'echo "Prerequisite file generation finished for trio: {offspring_id}"\n')
        f_sh.write(f'echo "End Date: $(date)"\n')

    # --- Submit the .sh to the server ---
    print(f"Submitting prerequisite job for trio {offspring_id}...")
    sub_cmd = f"sbatch {sh_script_name}"
    try:
        # Use subprocess.run to capture potential sbatch errors
        result = subprocess.run(sub_cmd, shell=True, check=True, capture_output=True, text=True)
        job_id = result.stdout.strip().split()[-1] # sbatch typically prints "Submitted batch job <ID>"
        print(f"Job submitted successfully. Job ID: {job_id}. Check logs in {trio_output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to submit job for {offspring_id}.")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")


################
# What you run #
################

if __name__ == "__main__":
    # Ensure the main phasing directory exists
    os.makedirs(direct_phasing, exist_ok=True)
    
    print("Starting SBATCH submission for prerequisite file generation pipeline for phasing.")
    for offspring_id, mother_id, father_id in trio_data.values():
        prepare_trio_phasing_files_sbatch(offspring_id=offspring_id, 
                                   mother_id=mother_id, 
                                   father_id=father_id, 
                                   denovo_dir=direct_denovo,
                                   phasing_dir=direct_phasing,
                                   ref_fasta_dir=ref_fasta_dir,
                                   vcf_files_dir=direct_vcf_files,
                                   window_flank=PHASING_WINDOW_FLANK,
                                   account_name=account, # From variable.py
                                   scratch_dir=scratch_dir, # From variable.py (though not used directly in this script)
                                   ref_genome_name=refGenome) # From variable.py
    print("\nAll prerequisite file generation jobs submitted.")
