## This script creates all the directories needed for the pipeline.
SPECIES=Odocoileus_virginianus
Path=/home/projects/rrg-shaferab/devan

############################################################################################################
mkdir $Path/$SPECIES/

mkdir $Path/$SPECIES/ref_fasta

mkdir $Path/$SPECIES/fastq_files
mkdir $Path/$SPECIES/fastq_files/catenate

mkdir $Path/$SPECIES/trimmed_seq

mkdir $Path/$SPECIES/bam_files
mkdir $Path/$SPECIES/bam_files/coverage
mkdir $Path/$SPECIES/bam_files/coverage/cov.log
mkdir $Path/$SPECIES/bam_files/inter_bam
mkdir $Path/$SPECIES/bam_files/map.log
mkdir $Path/$SPECIES/bam_files/merge.log
mkdir $Path/$SPECIES/bam_files/recal.log
mkdir $Path/$SPECIES/bam_files/summary
mkdir $Path/$SPECIES/bam_files/uniq_rmdup.log

mkdir $Path/$SPECIES/vcf_files
mkdir $Path/$SPECIES/vcf_files/inter_vcf
mkdir $Path/$SPECIES/vcf_files/back_com.log
mkdir $Path/$SPECIES/vcf_files/call.log
mkdir $Path/$SPECIES/vcf_files/combine.log
mkdir $Path/$SPECIES/vcf_files/gather.log

mkdir $Path/$SPECIES/vcf_handling
mkdir $Path/$SPECIES/vcf_handling/log_file
mkdir $Path/$SPECIES/vcf_handling/inter_vcf

mkdir $Path/$SPECIES/de_novo_mutation
mkdir $Path/$SPECIES/de_novo_mutation/log_file
mkdir $Path/$SPECIES/de_novo_mutation/inter_vcf
