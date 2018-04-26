#!/usr/bin/env python

#Remove sequences of fly lines that are complelety unusable (Ns at all positions)
#Nov 11th 2017
#Joo Hyun Im (ji72)


import sys

from collections import Counter

header = ''

for line in open(sys.argv[1]):
    if '>' in line:
        header = line
    elif line == '\n':
        pass
    elif line.count('N') < len(line)*0.1: #only pick lines that have less than 10% of the sequence as missing
        print header,
        print line,
