#!/bin/bash -x
#To trim these sequences and merge them:
#pushd .. && make clean && make && popd

#6 * 3 * 5 * 4
Z=( 26 )
M=( 0.02 )
N=( 0.87 )
X=( 0.125 )
M2=( 0.02 )
N2=( 0.9 )
Q=( 13)

rm trimmed_*

for z in "${Z[@]}"
do
    for m in "${M[@]}"
    do
	for n in "${N[@]}"
	do
	    for x in "${X[@]}"
	    do
		for q in "${Q[@]}"
		do
		    for m2 in "${M2[@]}"
		    do
			for n2 in "${N2[@]}"
			do
		    ../../SeqPrep \
			-f ./simSeq10k_1.fq \
			-r ./simSeq10k_2.fq \
			-1 ./simSeq_trimmed_1.fq.gz \
			-2 ./simSeq_trimmed_2.fq.gz \
			-Z ${z} \
			-M ${m} \
			-N ${n} \
			-X ${x} \
			-q ${q} \
			-n ${n2} \
			-m ${m2}
		    rm simSeq_trimmed_1.fq
		    rm simSeq_trimmed_2.fq
		    gunzip simSeq_trimmed_1.fq.gz
		    gunzip simSeq_trimmed_2.fq.gz
		    python ./simseq_trimmed_error_check.py simSeq10k_1.fq simSeq10k_2.fq simSeq_trimmed_1.fq simSeq_trimmed_2.fq > trimmed_M${m}_N${n}_X${x}_Z${z}_Q${q}_m${m2}_n${n2}.txt
			done
		    done
		done
	    done
	done
    done
done


#cp simSeq10k_1.fq simSeq10k_mcf_1.fq
#cp simSeq10k_1.fq simSeq10k_mcf_2.fq
#ea-utils-read-only/clipper/fastq-clipper simSeq10k_mcf_1.fq
#ea-utils-read-only/clipper/fastq-clipper -o simSeq10k_mcf_1.fq simSeq10k_1.fq AGATCGGAAGAGCGGTTCAG
#ea-utils-read-only/clipper/fastq-clipper -o simSeq10k_mcf_2.fq simSeq10k_2.fq AGATCGGAAGAGCGTCGTGT
#python ./simseq_trimmed_error_check.py simSeq10k_1.fq simSeq10k_2.fq simSeq10k_mcf_1.fq simSeq10k_mcf_2.fq >trimmed_M99_N99_X99_Z99_Q99_m99_n99.txt


python sens_vs_spec.py
#open result.html
