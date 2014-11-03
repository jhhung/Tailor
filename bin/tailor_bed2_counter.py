#! /usr/bin/env python
# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

'''
this script takes the BED2 format given by Tailor
and summarize the tailing events
'''
# bed2 format specification:
# f1-6: same as ordinary bed
# f7: original sequence as in input
# f8: sequence of the tail (* if no tail)
# f9: length of tail 
# f10: MAPQ string for mismatches (* if not mismatch)

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

