#!/bin/bash
##
## This script filters low quality reads
##
source ../../variable.py
##

##
###############################
## Construct the file clean.sh file:
less $path/$sp/raw_seq_dir.txt | while read a b c; do echo "SOAPnuke filter -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -1 $b -2 $c -G -Q 2 -l 10 -q 0.2 -E 60 -5 0 -M 2 -o $path/$ap/trimmed_seq/$(basename $b '_read_1.fq.gz') -C $(basename $b '.fq.gz').clean.fq.gz -D $(basename $c '.fq.gz').clean.fq.gz"; done > clean.sh
##
## Split the files to submit them:
bash submit_split.sh
##
chmod u+rwx clean_splitted.sh
bash clean_splitted.sh
##
## Now the clean files are created in $path/$sp/trimmed_seq
