## This script creates all the directories needed for the pipeline.
sp=Odocoileus_virginianus
path=/home/projects/rrg-shaferab/devan

############################################################################################################
mkdir $path/$sp/

mkdir $path/$sp/ref_fasta

mkdir $path/$sp/fastq_files
mkdir $path/$sp/fastq_files/catenate

mkdir $path/$sp/trimmed_seq

mkdir $path/$sp/bam_files
mkdir $path/$sp/bam_files/coverage
mkdir $path/$sp/bam_files/coverage/cov.log
mkdir $path/$sp/bam_files/inter_bam
mkdir $path/$sp/bam_files/map.log
mkdir $path/$sp/bam_files/merge.log
mkdir $path/$sp/bam_files/recal.log
mkdir $path/$sp/bam_files/summary
mkdir $path/$sp/bam_files/uniq_rmdup.log

mkdir $path/$sp/vcf_files
mkdir $path/$sp/vcf_files/inter_vcf
mkdir $path/$sp/vcf_files/back_com.log
mkdir $path/$sp/vcf_files/call.log
mkdir $path/$sp/vcf_files/combine.log
mkdir $path/$sp/vcf_files/gather.log

mkdir $path/$sp/vcf_handling
mkdir $path/$sp/vcf_handling/log_file
mkdir $path/$sp/vcf_handling/inter_vcf

mkdir $path/$sp/de_novo_mutation
mkdir $path/$sp/de_novo_mutation/log_file
mkdir $path/$sp/de_novo_mutation/inter_vcf
