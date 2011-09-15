#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
seqlens = defaultdict(int)

from sys import stdin

next_line_seq = False
count = 0
for line in stdin:
	count += 1
	if line.startswith("@"):
		count = 0
		next_line_seq = True
	if next_line_seq and count == 1:
		next_line_seq == False
		seqlens[len(line)] += 1
		
	

for (length,count) in sorted(seqlens.iteritems(), key=itemgetter(0),reverse=True):
	print("%d\t%d"%(length,count))
