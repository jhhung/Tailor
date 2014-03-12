#! /usr/bin/env python
'''
this script takes the BED2 format given by Tailor
and summarize the tailing events
'''
from __future__ import print_function
import sys
from collections import defaultdict

def main (filename):
	# counter counts all tail
	counter  = defaultdict (lambda: defaultdict (float) )
	# counter2 only count the no-tail or single nucleotide tails
	counter2 = defaultdict (lambda: defaultdict (float) )
	for line in open (filename):
		vector = line.strip().split()
		counter[int (vector[2])- int(vector[1])][vector[7]] += float(vector[3])/float(vector[4])
		if (len (vector[7]) == 1):
			counter2[int (vector[2])- int(vector[1])][vector[7]] += float(vector[3])/float(vector[4])
		else:
			counter2[int (vector[2])- int(vector[1])]["the_others"] += float(vector[3])/float(vector[4])
	for length in counter:
		for tails in counter[length]:
			print (length, tails, counter[length][tails], sep='\t', end='\n', file=sys.stdout)
	for length in counter2:
		for tails in counter2[length]:
			print (length, tails, counter2[length][tails], sep='\t', end='\n', file=sys.stderr)
	
if __name__ == '__main__' :
	if len (sys.argv) != 2:
		print ("usage:", sys.argv[0], "input.bed", sep=' ', end='\n', file=sys.stderr)
	else:
		main (sys.argv[1])

