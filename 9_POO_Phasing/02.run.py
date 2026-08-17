# -*- coding: utf-8 -*-
"""
This script runs the parent-of-origin phasing script (assign_parent_of_origin.py)
locally for each trio defined in the pedigree.
"""
#################
# What you need #
#################

# Packages:
import subprocess
import os
import pandas as pd
from variable import * # Import your variables.py (sp, path, refGenome)

# --- Directories ---
direct_phasing = "{}/{}/phasing/".format(path, sp) # Directory where prerequisite files are located
direct_bam_files = "{}/{}/bam_files/".format(path, sp) # Directory containing BAM files
ref_fasta_dir = "{}/{}/ref_fasta/".format(path, sp) # Directory for FASTA
assign_parent_script_path = os.path.join(path, sp, "phasing", "assign_parent_of_origin.py") # Path to the phasing script itself
direct = "{}/{}/".format(path, sp) # Base project directory (from variable.py context)

# --- Parameters for phasing (ensure consistency with prepare_phasing_prerequisites.py) ---
PHASING_WINDOW = 250 # Window size used for flanking regions in VCF subsetting

# --- Dictionary for trios ---
pedigree_path = '{}pedigree.ped'.format(direct) # Now 'direct' is explicitly defined
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


# --- The function to run the phasing job for a single trio locally ---
def run_phasing_job_locally(offspring_id, mother_id, father_id, 
                            phasing_base_dir, bam_files_dir, fasta_dir, 
                            window_size, ref_genome_name, assign_parent_script):
    """
    Runs the assign_parent_of_origin.py script locally for a single trio.
    """
    
    print(f"\n--- Running parent-of-origin phasing for trio: {offspring_id} ---")

    # Trio-specific directory (where prerequisites were generated and results will go)
    trio_output_dir = os.path.join(phasing_base_dir, offspring_id)
    os.makedirs(trio_output_dir, exist_ok=True) # Ensure it exists

    # --- Input file paths for assign_parent_of_origin.py ---
    vcf_file = os.path.join(trio_output_dir, f"dnms_flank_{offspring_id}.vcf.gz")
    dnm_sites_file = os.path.join(trio_output_dir, f"sites_{offspring_id}.txt")
    fasta_file = os.path.join(fasta_dir, f"{ref_genome_name}.fa")

    # --- Output file path for phasing results ---
    phasing_results_file = os.path.join(trio_output_dir, f"{offspring_id}_phasing_results.tsv")
    
    # Construct the command to run assign_parent_of_origin.py
    command = [
        "python", assign_parent_script,
        "--vcf", vcf_file,
        "--bam_dir", bam_files_dir,
        "--fasta", fasta_file,
        "--offspring", offspring_id,
        "--mother", mother_id,
        "--father", father_id,
        "--dnm_sites", dnm_sites_file,
        "--output", phasing_results_file,
        "--window", str(window_size)
    ]

    print(f"  Executing command for {offspring_id}:")
    print(f"    {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"  Stdout:\n{result.stdout.strip()}")
        if result.stderr:
            print(f"  Stderr:\n{result.stderr.strip()}")
        print(f"--- Phasing job completed for trio: {offspring_id} ---\n")
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Phasing job failed for trio {offspring_id}: {e}")
        print(f"  Command: {' '.join(command)}")
        print(f"  Stdout: {e.stdout}")
        print(f"  Stderr: {e.stderr}")
        print(f"  Stopping processing for trio {offspring_id} due to error.")
        # Do not return here if you want to proceed with other trios and aggregation
        # If you want to stop the entire script on first error, then 'exit(1)'
        # For aggregation, it's better to just skip this trio's data and continue.


################
# What you run #
################

if __name__ == "__main__":
    os.makedirs(direct_phasing, exist_ok=True)
    
    print("Starting local parent-of-origin phasing pipeline.")
    
    for offspring_id, mother_id, father_id in trio_data.values():
        run_phasing_job_locally(offspring_id=offspring_id, 
                                mother_id=mother_id, 
                                father_id=father_id, 
                                phasing_base_dir=direct_phasing,
                                bam_files_dir=direct_bam_files,
                                fasta_dir=ref_fasta_dir,
                                window_size=PHASING_WINDOW,
                                ref_genome_name=refGenome,
                                assign_parent_script=assign_parent_script_path)
        
    print("\nAll parent-of-origin phasing tasks completed.")

    # --- Aggregation and Ratio Calculation ---
    print("\nAggregating phasing results and calculating ratios...")
    summary_data = []

    for offspring_id, _, _ in trio_data.values():
        trio_results_file = os.path.join(direct_phasing, offspring_id, f"{offspring_id}_phasing_results.tsv")
        
        if os.path.exists(trio_results_file) and os.path.getsize(trio_results_file) > 0:
            try:
                # Read the individual phasing results TSV
                df_results = pd.read_csv(trio_results_file, sep='\t')
                
                # Using 'parent_of_origin' as the column name for parent assignments
                # And using 'mother'/'father' for comparisons
                if 'parent_of_origin' in df_results.columns:
                    maternal_dnm_count = (df_results['parent_of_origin'] == 'mother').sum()
                    paternal_dnm_count = (df_results['parent_of_origin'] == 'father').sum()
                    
                    # Calculate Paternal to Maternal ratio, handling division by zero
                    if maternal_dnm_count > 0: # If maternal DNMs are > 0, calculate ratio
                        paternal_maternal_ratio = paternal_dnm_count / maternal_dnm_count 
                    else: # If maternal DNMs are 0, ratio is effectively infinite
                        paternal_maternal_ratio = float('inf') 
                        
                    summary_data.append({
                        'Offspring_ID': offspring_id,
                        'Maternal_DNMs': maternal_dnm_count,
                        'Paternal_DNMs': paternal_dnm_count,
                        'Paternal_to_Maternal_Ratio': paternal_maternal_ratio # Changed column name
                    })
                else:
                    print(f"Warning: 'parent_of_origin' column not found in {trio_results_file}. Skipping ratio calculation for {offspring_id}.")
                    summary_data.append({
                        'Offspring_ID': offspring_id,
                        'Maternal_DNMs': 'N/A',
                        'Paternal_DNMs': 'N/A',
                        'Paternal_to_Maternal_Ratio': 'N/A'
                    })
            except pd.errors.EmptyDataError:
                print(f"Warning: {trio_results_file} is empty. Skipping ratio calculation for {offspring_id}.")
            except Exception as e:
                print(f"Error reading or processing {trio_results_file}: {e}. Skipping ratio calculation for {offspring_id}.")
        else:
            print(f"Warning: Phasing results file not found or empty for {offspring_id}: {trio_results_file}. Skipping ratio calculation.")

    if summary_data:
        final_summary_df = pd.DataFrame(summary_data)
        # Output to /home/devan/projects/rrg-shaferab/devan/Odocoileus_virginianus/phasing/
        final_output_path = os.path.join(direct_phasing, "parent_of_origin_summary_ratios.tsv")
        final_summary_df.to_csv(final_output_path, sep='\t', index=False)
        print(f"\nFinal summary ratios saved to: {final_output_path}")
    else:
        print("\nNo phasing results found or processed for summary calculation.")
