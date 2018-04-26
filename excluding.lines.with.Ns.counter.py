#!/usr/bin/env python
import sys

from collections import Counter

header = ''

for line in open(sys.argv[1]):
    if '>' in line:
        header = line
    elif line == '\n':
        pass
    elif line.count('N') < len(line)*0.5:
        print header,
        #print line.count('N'),
        #print len(line)*0.5,
        print line,
