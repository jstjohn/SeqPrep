#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
from string import rstrip
seqlens = defaultdict(int)

from sys import stdin
from sys import argv

sizes = dict()
outsizes = defaultdict(int)
i1 = open(argv[1],'r')
i2 = open(argv[2],'r')
o1 = open(argv[3],'r')
o2 = open(argv[4],'r')

for f in [i1,i2]:
	count = 0
	for line in f:
		if line.startswith("@") and (count == 4 or count == 0):
			count = 0
			ins_len = int(line.split('/')[1])
			sizes[rstrip(line)] = ins_len
		count += 1
i1.close()
i2.close()

for f in [o1,o2]:
	next_line_seq = False
	count = 0
	pid = ''
	for line in f:
		if line.startswith("@") and (count == 4 or count == 0):
			count = 0
			pid = rstrip(line)
			next_line_seq = True
		count += 1
	       	if next_line_seq and count == 2:
	       		next_line_seq == False
			outsizes[pid] = len(line)
o1.close()
o2.close()
		

FP = 0 # trimmed genomic bases
FN = 0 # left adapter behind
TN = 0 # left genomic bases behind
TP = 0 # trimmed adapter
min_len = 30
max_len = 100
#sens = tp/(tp+fn)
#spec = tn/(tn+fp)

for (pid,truelen) in sizes.iteritems():
	outlen = outsizes[pid]
	## if truelen >= min_len:
	if truelen < outlen:
		TP += max_len - outlen #trimmed difference between outlen and orilen
		FP += 0 #didn't trim any genomic bases
		TN += truelen #left length genomic bases
		FN += outlen - truelen #difference is adapter
	if outlen == truelen:
		TP += max_len - outlen
		FP += 0
		TN += truelen
		FN += 0 #got all adapter
	if truelen > outlen:
		TP += max_len - min(max_len, truelen) #trimmed all of adapter
		FP += min(max_len,truelen) - outlen #trimmed the difference of genomic
		TN += outlen #didnt trim outlen genomic bases
		FN += 0 #since the true length is > outlength there is no adapter left behind
	## else:
	## 	if outlen > 0:
	## 		FP += 0 # trimmed genomic bases
	## 		FN += max(outlen - truelen,0) # left adapter behind
	## 		TN += 0 # left genomic bases behind
	## 		TP += max_len - outlen # trimmed adapter
	## 	if outlen == 0:
	## 		FP += 0 # trimmed genomic bases
	## 		FN += 0 # left adapter behind
	## 		TN += 0 # left genomic bases behind
	## 		TP += max_len - truelen # trimmed adapter
		
			
			
print("Adapter Trimmed:%d" % (TP) )
print("Adapter Missed:%d" % (FN) )
print("Adapter Total:%d" % (FN+TP)  )
print("Genomic Bases Trimmed:%d" % (FP) )
print("Genomic Bases Total:%d" % (FP + TN) )
print("Adapter Trimming Sensitivity: %f" % (float(TP)/(float(TP)+float(FN))) )
print("Adapter Trimming Specificity: %f" % (float(TN)/(float(TN)+float(FP))) )
