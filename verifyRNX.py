#!/usr/bin/env python3

"""This script takes a RINEX file as input and verifies that it is of a
suitable length for inclusion in the NGCA. That is greater than 6 hrs but less
than 48 hrs.

Usage:

    verifyRNX <RINEX file>
"""

import re
import sys
import os
import datetime
import statistics


# Compile the regular expressions
p1 = re.compile(r'RINEX VERSION / TYPE')
p2_str = '^>\s\d{4}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}'
p2_str += '\s{1,2}\d{1,2}\.\d{,7}'
p2 = re.compile(p2_str)
p3_str = '^\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}'
p3_str += '\d{1,2}\s{1,2}\d{1,2}\.\d{,7}'
p3 = re.compile(p3_str)
p4 = re.compile('^\s\d\s')
p8 = re.compile('^\d{4}')

# Read in the RINEX file
rnx = sys.argv[1]
with open(rnx) as f:
    rnx_lines = f.readlines()

# Loop over the lines
epochs = []
deltas = []
for line in rnx_lines:
    line.strip()
    
    # Get the RINEX version
    if p1.search(line):
        version = float(line.split()[0])
        if version >= 3.0:
            rnx3 = True
        else:
            rnx3 = False

    # Turn the epoch lines into datetime objects and add to an array
    if rnx3:
        if p2.match(line):
            date_string = '{:4s} {:2s} {:2s} {:2s} {:2s} {:2s} {:6s}' \
                          .format(line[2:6], line[7:9], line[10:12],
                                  line[13:15],line[16:18], line[19:21],
                                  line[22:28])
            epoch = datetime.datetime.strptime(date_string, \
                    '%Y %m %d %H %M %S %f')
            epochs.append(epoch)
    elif not rnx3:
        if p3.match(line):
            date_string = '{:2s} {:2s} {:2s} {:2s} {:2s} {:2s} {:6s}' \
                          .format(line[1:3], line[4:6], line[7:9],line[10:12],
                                  line[13:15],line[16:18],line[19:25])
            if p4.match(date_string):
                parts = date_string.split()
                for i in range(len(parts)):
                    if len(parts[i]) == 1:
                        parts[i] = '0' + parts[i]
                date_string = ' '.join(parts)
            if date_string[15:17] == '60':
                tmp = date_string[:15] + '00' + date_string[17:]
                date_string = tmp
            if date_string[18:25] != '000000':
                seconds = float(date_string[15:17] + '.' + date_string[18:25])
            epoch = datetime.datetime.strptime(date_string, \
                    '%y %m %d %H %M %S %f')
            epochs.append(epoch)
        
# Check that the DOY in file name is right
doy = epochs[0].strftime('%j')
if rnx[4:7] != doy:
    sys.exit('Incorrect DOY in RINEX file name')

# Calculate the sampling interval
for i in range(1, len(epochs)):
    delta_t = epochs[i] - epochs[i-1]
    deltas.append(delta_t.seconds)
sampling = statistics.median(deltas)
        
# Calculate the duration of the RINEX file
duration = len(deltas) * sampling

# 05:59:00 hours in seconds
if duration < 21540:
    sys.exit('RINEX duration is less than 6 hours (' + str(duration) + ')')

# 48 hours in seconds    
if duration > 172800:
    sys.exit('RINEX duration is greater than 48 hours (' + str(duration) + ')')

# Rename file
os.rename(rnx, rnx.upper())
