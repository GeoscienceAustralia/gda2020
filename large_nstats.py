#!/usr/bin/env python3

"""A script to find the largest nstats in a DynAdjust .adj file 
"""

import re
import sys


# Create the regular expression
p1 = re.compile('^Adjusted Measurements')
p2 = re.compile('^Adjusted Coordinates')

# Open the input file and read in the lines
if len(sys.argv) != 2:
    print('Please specify a DynAdjust .adj file') 
    sys.exit()
infile = sys.argv[1]
f = open(infile)
lines = f.readlines()

# Create list to hold the adjusted measurement data
adj_msrs = []
header = False
go = False
test = False
n = 0
for line in lines:
    line = line.rstrip()
    if p1.match(line):
        header = True
    if p2.match(line):
        break
    if header:
        if n == 4:
            header = False
            go = True
        if n < 4:
            n += 1
            continue
    elif go:
        if line != '':
            adj_msrs.append(line)

# Loop over the adj_msrs create a list of nstats 
all_nstats = []
for adj_msr in adj_msrs:
    nstat = adj_msr[157:168]
    nstat = nstat.strip()
    if nstat != '' and nstat != '*':
        all_nstats.append(float(nstat))

# Sort the list and create a second list with the 10 biggest positive nstats
# and the 10 biggest negatve nstats 
all_nstats.sort()
big_pos_nstats = all_nstats[-10:] 
big_neg_nstats = all_nstats[0:10]
big_nstats = big_pos_nstats + big_neg_nstats
abs_nstats = list(map(abs, big_nstats))
all_values = list(zip(big_nstats, abs_nstats))
all_values.sort(key=lambda a: a[1])
values = all_values[-10:]
nstats, abs_nstats = zip(*values)
print('The biggest nstats are:')
for nstat in nstats:
    print(nstat)
