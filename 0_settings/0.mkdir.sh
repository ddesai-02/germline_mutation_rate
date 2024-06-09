## This script creates all the directories needed for the pipeline.
sp=Odocoileus_virginianus
Path=/home/projects/rrg-shaferab/devan

############################################################################################################
mkdir $Path/$sp/

mkdir $Path/$sp/ref_fasta

mkdir $Path/$sp/fastq_files
mkdir $Path/$sp/fastq_files/catenate

mkdir $Path/$sp/trimmed_seq

mkdir $Path/$sp/bam_files
mkdir $Path/$sp/bam_files/coverage
mkdir $Path/$sp/bam_files/coverage/cov.log
mkdir $Path/$sp/bam_files/inter_bam
mkdir $Path/$sp/bam_files/map.log
mkdir $Path/$sp/bam_files/merge.log
mkdir $Path/$sp/bam_files/recal.log
mkdir $Path/$sp/bam_files/summary
mkdir $Path/$sp/bam_files/uniq_rmdup.log

mkdir $Path/$sp/vcf_files
mkdir $Path/$sp/vcf_files/inter_vcf
mkdir $Path/$sp/vcf_files/back_com.log
mkdir $Path/$sp/vcf_files/call.log
mkdir $Path/$sp/vcf_files/combine.log
mkdir $Path/$sp/vcf_files/gather.log

mkdir $Path/$sp/vcf_handling
mkdir $Path/$sp/vcf_handling/log_file
mkdir $Path/$sp/vcf_handling/inter_vcf

mkdir $Path/$sp/de_novo_mutation
mkdir $Path/$sp/de_novo_mutation/log_file
mkdir $Path/$sp/de_novo_mutation/inter_vcf
