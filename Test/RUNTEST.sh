#!/bin/bash -x
#To trim these sequences and merge them:
#pushd .. && make clean && make && popd
../SeqPrep -6 \
-f multiplex_bad_contam_1.fq.gz \
-r multiplex_bad_contam_2.fq.gz \
-A GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-B CAAGCAGAAGACGGCATACGAGA \
-1 pe_bad_contam_trimmed_1.fastq.gz \
-2 pe_bad_contam_trimmed_2.fastq.gz \
-s pe_bad_contam_merged.fastq.gz \
-E alignments.txt.gz

