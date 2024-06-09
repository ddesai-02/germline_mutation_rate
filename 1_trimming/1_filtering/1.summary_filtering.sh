#!/bin/bash
##
## This script runs the summary_SOP for all sequences
##
source ../variable.py
##

##
## Remove the output of cleanning:
rm submit_clean*.out
## Again build a bash with all the bash for that
ls -d $path/$sp/trimmed_seq/*/ | while read p; do echo "bash summary_SOAPnuke.sh $path/$sp/trimmed_seq/$(basename "$p" '/') $path/$sp/trimmed_seq"; done > run_summary_SOAPnuke.sh
##
##
chmod u+rwx run_summary_SOAPnuke.sh
bash run_summary_SOAPnuke.sh
##
## Now you should have summary_SOAPnuke.txt in the output directory : $PATH/$SP/trimmed_seq/
##
## Write the directories of the trimmed files
less $path/$sp/raw_seq_dir.txt | while read a b c; do echo -e "$a\t$(echo $b | sed -r "s:.fq.gz:.clean.fq.gz:g; s:fastq_files:trimmed_seq/$(basename "$b" '_read_1.fq.gz'):g;")\t$(echo $c | sed -r "s:.fq.gz:.clean.fq.gz:g; s:fastq_files:trimmed_seq/$(basename "$b" '_read_1.fq.gz'):g;")";done > $path/$sp/clean_seq_dir.txt
