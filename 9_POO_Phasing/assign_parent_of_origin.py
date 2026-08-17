#!/usr/bin/env python3

import pysam
from cyvcf2 import VCF
import argparse
from collections import defaultdict

# Function to parse and validate command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign parent-of-origin to DNMs using informative SNPs and read-backed phasing"
    )
    parser.add_argument("--vcf", required=True,
                        help="Trio VCF with DNMs ±50kb")
    parser.add_argument("--bam_dir", required=True,
                        help="Directory containing BAM files")
    parser.add_argument("--fasta", required=True,
                        help="Path to reference FASTA")
    parser.add_argument("--offspring", required=True,
                        help="Sample name for offspring in VCF and BAM naming")
    parser.add_argument("--mother", required=True,
                        help="Sample name for mother in VCF")
    parser.add_argument("--father", required=True,
                        help="Sample name for father in VCF")
    parser.add_argument("--window", type=int, default=50000,
                        help="Window size to search for informative SNPs around each DNM")
    parser.add_argument("--dnm_sites", required=True, # Made required for this iteration strategy
                        help="TSV file listing specific DNM sites to analyze")
    parser.add_argument("--output", required=False,
                        help="Path to output TSV file for parent-of-origin assignments")
    parser.add_argument("--debug_dnm", required=False,
                        help="Chrom:Pos of a DNM to debug read collection (e.g., 'chr1:12345')")
    return parser.parse_args()

# Function: Convert a genotype string to a set of its constituent nucleotide bases,
# supporting multi-allelic sites.
def gt_to_bases(gt, ref, alt_list): # alt_list is now a list of ALT alleles
    """
    Convert a genotype string to a set of possible nucleotide bases for multi-allelic sites.
    gt: genotype string (e.g., "0/0", "0/1", "0/2", "1/2")
    ref: Reference allele (string)
    alt_list: List of alternate alleles (e.g., ['T', 'G'])
    """
    bases = set()
    alleles_indices = gt.split('/')
    
    # Create a mapping from genotype index to actual base
    allele_map = {0: ref}
    for i, a in enumerate(alt_list):
        allele_map[i + 1] = a # 1 corresponds to first ALT, 2 to second, etc.

    for allele_idx_str in alleles_indices:
        try:
            allele_idx = int(allele_idx_str)
            if allele_idx in allele_map:
                bases.add(allele_map[allele_idx])
            else: # Handle cases where allele_idx is out of bounds for alt_list
                pass 
        except ValueError: # Handles "./." or other non-integer genotypes
            pass # Skip invalid alleles
    return bases

# Function: Determine if a trio SNP is potentially informative for phasing.
def is_informative(mom_gt, dad_gt, child_gt, debug_dnm_key=None, current_snp_pos=None):
    """
    Determine if a trio SNP is potentially informative for phasing.
    Simplified: child must be heterozygous, and parents must have genotypes (not missing).
    debug_dnm_key: The 'chrom:pos' string of the DNM being debugged, or None.
    current_snp_pos: The 1-based genomic position of the current SNP being evaluated.
    """
    if debug_dnm_key:
        print(f"  [DEBUG] is_informative called for SNP {current_snp_pos}: Mom GT={mom_gt}, Dad GT={dad_gt}, Child GT={child_gt}")

    if "." in (mom_gt, dad_gt, child_gt):
        if debug_dnm_key:
            print(f"    [DEBUG] SNP {current_snp_pos}: Skipping due to missing genotype.")
        return False
    
    # Child must be heterozygous at the informative SNP position to be able to phase a DNM in them.
    if child_gt not in {"0/1", "1/0", "0/2", "1/2", "2/1", "1/3", "3/1"}: # Expanded to cover multi-allelic child het
        if debug_dnm_key:
            print(f"    [DEBUG] SNP {current_snp_pos}: Skipping due to child not heterozygous ({child_gt}).")
        return False
    
    if debug_dnm_key:
        print(f"    [DEBUG] SNP {current_snp_pos}: Potentially informative (child is het).")
    return True # If child is het and parents have genotypes, it's potentially informative

# Function: Retrieve genotype (e.g., "0/1") using precomputed sample index
def get_genotype(record, sample_idx):
    """
    Retrieves the genotype string for a given sample index from a VCF record.
    record: A cyvcf2 VCF record object.
    sample_idx: The integer index of the sample in the VCF.
    Returns a genotype string (e.g., "0/1", "1/1", "0/0", or "./." for missing).
    """
    gt = record.genotypes[sample_idx]
    if len(gt) >= 2:
        return f"{gt[0]}/{gt[1]}"
    else:
        return "./."

# Function: Fetch reads overlapping given variant sites and record bases seen per-read
def fetch_reads_with_sites(bam_path, chrom, sites, debug_dnm_key=None):
    """
    Fetches reads from a BAM file that overlap specified sites and records the base
    observed at each site for each read.
    bam_path: Path to the BAM file.
    chrom: Chromosome name.
    sites: A list of 1-based genomic positions (integers) to fetch reads for.
    debug_dnm_key: The 'chrom:pos' string of the DNM being debugged, or None.
    Returns a defaultdict where keys are read names and values are sets of (site, base) tuples.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    reads_with_sites = defaultdict(set)
    
    if debug_dnm_key:
        print(f"  [DEBUG] fetch_reads_with_sites: Requesting reads for sites {sites} on {chrom}")

    # Process sites in chunks or individual fetches for robust debugging
    for site in sites:
        # Pysam coordinates are 0-based, VCF is 1-based, so site - 1 for fetch
        try:
            for read in bam.fetch(chrom, site - 1, site):
                if not read.is_unmapped and not read.is_duplicate and read.query_sequence is not None:
                    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                        if rpos == site - 1: # Check if the reference position matches our site (0-based)
                            reads_with_sites[read.query_name].add((site, read.query_sequence[qpos]))
        except ValueError as e:
            if "coordinate out of range" in str(e):
                if debug_dnm_key:
                    print(f"    [DEBUG] Skipping fetch for {chrom}:{site} as coordinates are out of range.")
            else:
                raise e # Re-raise other unexpected ValueErrors
    
    if debug_dnm_key:
        print(f"  [DEBUG] fetch_reads_with_sites: Collected {len(reads_with_sites)} unique reads covering requested sites.")
        for read_id, sites_covered in reads_with_sites.items():
            print(f"    [DEBUG] Read {read_id} covers: {sites_covered}")

    bam.close()
    return reads_with_sites

# Function: Compare reads covering both DNM and informative SNP to count parental support
def assign_parent(dnm_site, dnm_base, inf_site, mom_gt_snp, dad_gt_snp, ref_base_snp, alt_list_snp, read_info, debug_dnm_key=None):
    """
    Compares reads covering both the DNM and an informative SNP to tally parental support.
    dnm_site: 1-based position of the DNM.
    dnm_base: The alternate allele of the DNM that is being phased.
    inf_site: 1-based position of the informative SNP.
    mom_gt_snp: Genotype string of the mother at the informative SNP (e.g., "0/1").
    dad_gt_snp: Genotype string of the father at the informative SNP (e.g., "0/0").
    ref_base_snp: Reference allele at the informative SNP.
    alt_list_snp: List of alternate alleles at the informative SNP (e.g., ['T', 'G']).
    read_info: Dictionary of reads covering sites, as returned by fetch_reads_with_sites.
    debug_dnm_key: The 'chrom:pos' string of the DNM being debugged, or None.
    Returns a tuple (mom_support, dad_support) indicating counts of reads supporting each parent.
    """
    mom_support = 0
    dad_support = 0

    mom_alleles_at_snp = gt_to_bases(mom_gt_snp, ref_base_snp, alt_list_snp)
    dad_alleles_at_snp = gt_to_bases(dad_gt_snp, ref_base_snp, alt_list_snp)

    if debug_dnm_key:
        print(f"    [DEBUG] assign_parent for inf_SNP {inf_site}: Mom alleles={mom_alleles_at_snp}, Dad alleles={dad_alleles_at_snp}")

    for read_id, sites_covered_by_read in read_info.items():
        site_bases = dict(sites_covered_by_read) 
        
        if dnm_site in site_bases and inf_site in site_bases:
            read_dnm_base = site_bases[dnm_site]
            read_inf_base = site_bases[inf_site]

            if debug_dnm_key:
                print(f"      [DEBUG] Read {read_id}: DNM base={read_dnm_base} (expected {dnm_base}), Inf SNP base={read_inf_base}")

            if read_dnm_base == dnm_base: 
                # Prioritized logic for assigning parental support
                # 1. Check for alleles *unique* to the mother (triple-het, or other unique allele cases)
                if read_inf_base in mom_alleles_at_snp and read_inf_base not in dad_alleles_at_snp:
                    mom_support += 1
                    if debug_dnm_key:
                        print(f"        [DEBUG] Read {read_id}: MATERNAL unique allele support. Mom_support={mom_support}")
                # 2. Check for alleles *unique* to the father (triple-het, or other unique allele cases)
                elif read_inf_base in dad_alleles_at_snp and read_inf_base not in mom_alleles_at_snp:
                    dad_support += 1
                    if debug_dnm_key:
                        print(f"        [DEBUG] Read {read_id}: PATERNAL unique allele support. Dad_support={dad_support}")
                # 3. Check for shared alleles where one parent is homozygous and the other is heterozygous.
                elif (len(mom_alleles_at_snp) == 1 and read_inf_base in mom_alleles_at_snp and len(dad_alleles_at_snp) > 1 and read_inf_base in dad_alleles_at_snp):
                    mom_support += 1
                    if debug_dnm_key:
                        print(f"        [DEBUG] Read {read_id}: MATERNAL homozygous/heterozygous support (mom homozygous). Mom_support={mom_support}")
                elif (len(dad_alleles_at_snp) == 1 and read_inf_base in dad_alleles_at_snp and len(mom_alleles_at_snp) > 1 and read_inf_base in mom_alleles_at_snp):
                    dad_support += 1
                    if debug_dnm_key:
                        print(f"        [DEBUG] Read {read_id}: PATERNAL homozygous/heterozygous support (dad homozygous). Dad_support={dad_support}")
                else:
                    if debug_dnm_key:
                        print(f"        [DEBUG] Read {read_id}: No clear parental support from this read/SNP combination.")
        else:
            if debug_dnm_key:
                print(f"      [DEBUG] Read {read_id}: Does not cover both DNM ({dnm_site}) and Inf SNP ({inf_site}).")
                        
    return mom_support, dad_support

def main():
    args = parse_args()
    vcf = VCF(args.vcf)

    # Map sample names to VCF columns
    try:
        mother_idx = vcf.samples.index(args.mother)
        father_idx = vcf.samples.index(args.father)
        offspring_idx = vcf.samples.index(args.offspring)
    except ValueError as e:
        print(f"Error: Sample name not found in VCF header: {e}")
        return

    # Load all DNM sites from the --dnm_sites file
    dnm_sites_to_process = []
    if args.dnm_sites: # This argument is now required
        with open(args.dnm_sites) as f:
            for line in f:
                chrom, pos_str = line.strip().split("\t")
                dnm_sites_to_process.append((chrom, int(pos_str)))
        print(f"Loaded {len(dnm_sites_to_process)} DNM sites from --dnm_sites file.")

    output_lines = []
    bam_path = f"{args.bam_dir}/{args.offspring}_sorted.merged.addg.uniq.rmdup.bam"

    # Try to open BAM file once to check for errors
    try:
        temp_bam = pysam.AlignmentFile(bam_path, "rb")
        temp_bam.close()
    except FileNotFoundError:
        print(f"Error: BAM file not found at {bam_path}. Please check --bam_dir and --offspring.")
        return
    except Exception as e:
        print(f"Error opening BAM file {bam_path}: {e}")
        return

    # Main processing loop: Iterate over the specified DNM sites and query the VCF for each
    for dnm_chrom, dnm_pos in dnm_sites_to_process:
        dnm_key = f"{dnm_chrom}:{dnm_pos}"

        # If --debug_dnm is specified, skip all other DNMs
        debug_dnm_active = (args.debug_dnm and args.debug_dnm == dnm_key)
        if args.debug_dnm and args.debug_dnm != dnm_key:
            continue

        print(f"[DEBUG] Attempting to process DNM: {dnm_key}")

        # Explicitly query the VCF for this specific DNM site using the VCF object as a function
        # This returns an iterator, so convert to list to check if records exist
        recs_at_dnm_pos = list(vcf(f"{dnm_chrom}:{dnm_pos}-{dnm_pos}")) # Corrected to vcf(region_string)

        if not recs_at_dnm_pos:
            print(f"[DEBUG] DNM {dnm_key} from dnm_sites.txt not found in VCF (no record at this exact position). Skipping.")
            output_lines.append(f"{dnm_chrom}\t{dnm_pos}\tnot_found_in_vcf")
            continue
        
        # We expect exactly one record for a simple DNM at a site
        # If more than one, pick the first one, or handle complex cases if necessary.
        # For typical SNV DNMs, there should be only one.
        rec = recs_at_dnm_pos[0]

        dnm_ref = rec.REF
        dnm_alt_list = rec.ALT # This is a list of ALT alleles from the VCF record

        # Filter 1: Only consider SNVs that are potentially de novo (e.g., allele count of 1 in the cohort)
        # Ensure it's a single nucleotide variant (SNV) and not an indel for this logic.
        # Note: rec.ALT is a list, so check len(rec.ALT) == 1 for single ALT allele SNV.
        if not rec.is_snp or len(rec.ALT) != 1 or rec.INFO.get("AC") != 1:
            print(f"[DEBUG] Skipping DNM {dnm_key}: Not a single-allele SNP, multi-allelic, or AC != 1. (is_snp={rec.is_snp}, len(ALT)={len(rec.ALT)}, AC={rec.INFO.get('AC')})")
            output_lines.append(f"{dnm_chrom}\t{dnm_pos}\tfiltered_out_by_variant_type")
            continue
        
        # Now we know it's a single ALT, so dnm_alt is safe.
        dnm_alt = dnm_alt_list[0]

        # Filter 2: A DNM should be heterozygous in the child (0/1 or 1/0, or multi-allelic hets like 0/2)
        child_gt = get_genotype(rec, offspring_idx)
        # Using a more robust check for heterozygosity, assuming 0/1, 1/0, 0/2, 2/0, 1/2, 2/1 etc.
        # This checks if the two alleles are different, and neither is missing.
        # It handles multi-allelic sites where child might be 0/2, 1/2 etc.
        # For DNMs, we expect 0/1 or 1/0 for single-ALT.
        if child_gt not in {"0/1", "1/0"}:
            print(f"[DEBUG] Skipping DNM {dnm_key}: Child genotype is not heterozygous (0/1 or 1/0): {child_gt}.")
            output_lines.append(f"{dnm_chrom}\t{dnm_pos}\tfiltered_out_by_child_gt")
            continue

        print(f"Analyzing DNM at {dnm_key} {dnm_ref}>{dnm_alt}")

        # Gather potentially informative SNPs around DNM
        nearby_snps = []
        try:
            # Query VCF for variants within the window around the DNM
            # Ensure query range is valid (pos-args.window might be < 1)
            query_start = max(1, dnm_pos - args.window) # Use dnm_pos here
            query_end = dnm_pos + args.window # Use dnm_pos here
            for snp in vcf(f"{dnm_chrom}:{query_start}-{query_end}"): # Corrected to vcf(region_string)
                # Skip if it's not a SNP or if it's the DNM site itself
                # Also skip multi-allelic SNPs for informative sites to keep logic simpler
                if not snp.is_snp or snp.POS == dnm_pos or len(snp.ALT) != 1:
                    if debug_dnm_active: # Only print this verbose debug if debug_dnm is active
                        if not snp.is_snp:
                            print(f"  [DEBUG] Skipping nearby variant {snp.CHROM}:{snp.POS}: Not a SNP.")
                        elif snp.POS == dnm_pos:
                            print(f"  [DEBUG] Skipping nearby variant {snp.CHROM}:{snp.POS}: It's the DNM site itself.")
                        elif len(snp.ALT) != 1:
                            print(f"  [DEBUG] Skipping nearby variant {snp.CHROM}:{snp.POS}: Multi-allelic informative SNP.")
                    continue
                
                mom_gt_snp = get_genotype(snp, mother_idx)
                dad_gt_snp = get_genotype(snp, father_idx)
                child_gt_snp = get_genotype(snp, offspring_idx)
                
                # Check if the SNP is informative for phasing (child is het, parents have genotypes)
                # Pass snp.ALT (the list) to is_informative if it needs multi-allelic context (not for this version)
                if is_informative(mom_gt_snp, dad_gt_snp, child_gt_snp, dnm_key if debug_dnm_active else None, snp.POS):
                    # Store SNP position, REF, the *list* of ALTs, and parental genotypes
                    nearby_snps.append((snp.POS, snp.REF, snp.ALT, mom_gt_snp, dad_gt_snp)) # snp.ALT is a list here!
        except Exception as e:
            if debug_dnm_active: # Only print this verbose debug if debug_dnm is active
                print(f"[DEBUG] Error querying VCF for nearby SNPs around {dnm_key}: {e}")
            pass # Continue processing even if there's an error fetching nearby SNPs

        if not nearby_snps:
            print(f"No *potentially* informative SNPs within ±{args.window}bp of {dnm_key}, undetermined.")
            output_lines.append(f"{dnm_chrom}\t{dnm_pos}\tundetermined_no_informative_snps")
            continue

        if debug_dnm_active: # Only print this verbose debug if debug_dnm is active
            print(f"  [DEBUG] Found {len(nearby_snps)} potentially informative SNPs around {dnm_key}:")
            for s_pos, s_ref, s_alt_list, s_mom_gt, s_dad_gt in nearby_snps:
                print(f"    [DEBUG] SNP: {dnm_chrom}:{s_pos} {s_ref}>{s_alt_list} MomGT:{s_mom_gt} DadGT:{s_dad_gt}")

        # Sort anchors by proximity to the DNM (closer SNPs might be more reliable)
        nearby_snps.sort(key=lambda x: abs(x[0]-dnm_pos))

        # Read collection: Fetch reads overlapping the DNM and all candidate informative SNPs
        sites_to_fetch = [dnm_pos] + [s[0] for s in nearby_snps]
        read_info = fetch_reads_with_sites(bam_path, dnm_chrom, sites_to_fetch, dnm_key if debug_dnm_active else None)

        # Debug output if a specific DNM is requested for debugging
        if debug_dnm_active: # Only print this verbose debug if debug_dnm is active
            print(f"Debug {dnm_key}: Read information collected: {read_info}")

        # Tally support for maternal and paternal origin
        total_mom = total_dad = 0
        # When iterating through nearby_snps, we now get 'alt_list_snp' which is a list
        for snp_pos, ref_base_snp, alt_list_snp, mom_gt_snp, dad_gt_snp in nearby_snps:
            # Pass the list of ALT alleles to assign_parent
            m, d = assign_parent(dnm_pos, dnm_alt, snp_pos, mom_gt_snp, dad_gt_snp, ref_base_snp, alt_list_snp, read_info, dnm_key if debug_dnm_active else None)
            total_mom += m
            total_dad += d

        # Call origin based on support thresholds
        if total_mom >= 10 and total_dad <= 2:
            call = 'mother'
        elif total_dad >= 10 and total_mom <= 2:
            call = 'father'
        else:
            call = 'undetermined'
        print(f"Total mom: {total_mom}, dad: {total_dad} -> {call}")
        output_lines.append(f"{dnm_chrom}\t{dnm_pos}\t{call}")

    # Write output to file if specified
    if args.output:
        with open(args.output,'w') as out:
            out.write("chrom\tpos\tparent_of_origin\n")
            for l in output_lines:
                out.write(l+"\n")

if __name__ == '__main__':
    main()
