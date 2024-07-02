# Enter species name and reference genome name without the ".fa" extension:
sp="Odocoileus_virginianus"                                 # Species name
refGenome="odocoileus_virginianus_NEW"                      # Name of reference genome
path="/home/devan/projects/rrg-shaferab/devan"              # Your project directory
nb_scaff=482                                                # Number of scaffolds
account="rrg-shaferab"                                      # Your SBATCH account
scratch_dir="/home/devan/scratch/"                          # Your scratch directory
ntfy_server="http://ntfy.metri.cc/germline_mutation_rate"   #ntfy server, set to empty if you're not using it

## ONLY FOR 6_DMNs
GQ_lim="60"
DP_min="0.5*mean_dp"
DP_max="2*mean_dp"
AB_max="0.7"
AB_min="0.3"
