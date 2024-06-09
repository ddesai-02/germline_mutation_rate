#!/usr/bin/env bash
# The reference genome should be indexed with:
source ../variable.py

bwa index -p $refGenome -a bwtsw $refGenome.fa
java -jar picard.jar CreateSequenceDictionary R=$refGenome.fa O=$refGenome.dict
samtools faidx $refGenome.fa
