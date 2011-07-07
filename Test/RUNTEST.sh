#!/bin/bash -x
#To trim these sequences and merge them:
#pushd .. && make clean && make && popd
../SeqPrep -6 \
-f ./data/multiplex_bad_contam_1.fq.gz \
-r ./data/multiplex_bad_contam_2.fq.gz \
-A GATCGGAAGAGCACACGTCT \
-B AGATCGGAAGAGCGTCGT \
-1 ./out/pe_bad_contam_merged_1.fastq.gz \
-2 ./out/pe_bad_contam_merged_2.fastq.gz \
-s ./out/pe_bad_contam_merged_s.fastq.gz \
-E ./info/alignments_merged.txt.gz

../SeqPrep -6 \
-f ./data/multiplex_bad_contam_1.fq.gz \
-r ./data/multiplex_bad_contam_2.fq.gz \
-A GATCGGAAGAGCACACGTCT \
-B AGATCGGAAGAGCGTCGT \
-1 ./out/pe_bad_contam_trimmed_1.fastq.gz \
-2 ./out/pe_bad_contam_trimmed_2.fastq.gz \
-E ./info/alignments_trimmed.txt.gz

prog=gzcat
$prog ./out/pe_bad_contam_trimmed_1.fastq.gz | python seqlens.py > ./info/pe_bad_contam_trimmed_1.lenhist.txt
$prog ./out/pe_bad_contam_trimmed_2.fastq.gz | python seqlens.py > ./info/pe_bad_contam_trimmed_2.lenhist.txt
$prog ./out/pe_bad_contam_merged_1.fastq.gz | python seqlens.py > ./info/pe_bad_contam_merged_1.lenhist.txt
$prog ./out/pe_bad_contam_merged_2.fastq.gz | python seqlens.py > ./info/pe_bad_contam_merged_2.lenhist.txt
$prog ./out/pe_bad_contam_merged_s.fastq.gz | python seqlens.py > ./info/pe_bad_contam_merged_s.lenhist.txt
