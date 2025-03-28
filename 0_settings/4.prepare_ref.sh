#!/usr/bin/env bash
# The reference genome should be indexed with:
source ../variable.py

echo $path/$sp/ref_fasta/$refGenome
bwa index -p $path/$sp/ref_fasta/$refGenome -a bwtsw $path/$sp/ref_fasta/$refGenome.fa
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$path/$sp/ref_fasta/$refGenome.fa O=$path/$sp/ref_fasta/$refGenome.dict
samtools faidx $path/$sp/ref_fasta/$refGenome.fa
