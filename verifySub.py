#!/usr/bin/env python3

"""This script verifies that the RINEX files in the current directory are
suitable for inclusion in the NGCA.

It runs the following checks:
    1) Antenna type is supported and height is a float 
        - this can be removed when turned into a function for the S3 bucket    
    2) Observation length (greater than 6 hrs but less than 48 hrs); 
    3) Sample interval (should be 30 seconds);
        - this can be removed when turned into a function for the S3 bucket    
    4) The file is from 1 June, 1994 or later no IGS products before that time;
        - this can be removed when turned into a function for the S3 bucket    
    5) DOY, i.e., that the DOY in the filename matches the first DOY of the data

RinexAntLs.txt is created from the information contained in the RINEX
headers and the log file lists all the RINEX files that were deleted
"""

import sys

import logging
import re
import os
import datetime
import statistics


# Set up logging
logging.basicConfig(
    filename = 'verifySub.log',
    level = logging.INFO
)

# Read the supported antennas into a set
antTypes = set()
with open('/home/fedora/antTypes.dat', 'r') as f:
    for line in f:
        cols = line.split()
        antTypes.add(cols[0])
        
# Compile the regular expressions
p1 = re.compile('^\w{8}\.\d{2}o$', re.I)
p2 = re.compile(r'RINEX VERSION / TYPE')
p3 = re.compile(r'ANT # / TYPE')
p4 = re.compile(r'ANTENNA: DELTA H/E/N')
p5 = re.compile('^\d+\.?\d+$')
p6a_str = '^>\s\d{4}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}'
p6a_str += '\s{1,2}\d{1,2}\.\d{,7}'
p6a = re.compile(p6a_str)
p6b_str = '^\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}\d{1,2}\s{1,2}'
p6b_str += '\d{1,2}\s{1,2}\d{1,2}\.\d{,7}'
p6b = re.compile(p6b_str)
p7 = re.compile('^\s\d\s')
p8 = re.compile('^\d{4}')

# Create a datetime object for the earliest IGS products
igs_start = datetime.datetime(1994, 6, 1, 0, 0, 0)

# Loop over the RINEX files
fout1 = open('RinexAntLs.txt', 'w')
fout2 = open('durations.dat', 'w')
stns = set()
num_rnx = 0
files = os.listdir('.')
for f in sorted(files):
    if p1.match(f):
        epochs = []
        deltas = []
        version = False
        antenna = False
        height = False
        good = True
        for line in open(f, 'r'):
            line.strip()

            # Get the RINEX version
            if version == False:
                if p2.search(line):
                    version = True
                    version = float(line.split()[0])
                    if version >= 3.0:
                        rnx3 = True
                    else:
                        rnx3 = False
            
            # Check that the antenna is supported
            if antenna == False:
                if p3.search(line):
                    antenna = True
                    antDome = line[20:40].strip()
                    ant = antDome.split()[0]
                    if ant not in antTypes:
                        logging.info(f + ': unsupported antenna type')
                        good = False
                        break
            
            # Check that the height is a number
            if height == False:
                if p4.search(line):
                    height = True
                    height = line[0:14].strip()
                    if not p5.match(height):
                        logging.info(f + ': antenna heights must be a number')
                        good = False
                        break
            
            # Turn the epoch lines into datetime objects and add to an array
            if rnx3:
                if p6a.match(line):
                    date_string = '{:4s} {:2s} {:2s} {:2s} {:2s} {:2s} {:6s}' \
                            .format(line[2:6], line[7:9], line[10:12],
                                    line[13:15], line[16:18], line[19:21],
                                    line[22:28])
                    epoch = datetime.datetime.strptime(date_string, \
                            '%Y %m %d %H %M %S %f')
                    epochs.append(epoch)
            elif not rnx3:
                if p6b.match(line):
                    date_string = '{:2s} {:2s} {:2s} {:2s} {:2s} {:2s} {:6s}' \
                            .format(line[1:3], line[4:6], line[7:9],
                                    line[10:12], line[13:15], line[16:18],
                                    line[19:25])
                    if p7.match(date_string):
                        parts = date_string.split()
                        for i in range(len(parts)):
                            if len(parts[i]) == 1:
                                parts[i] = '0' + parts[i]
                        date_string = ' '.join(parts)
                    if date_string[15:17] == '60':
                        tmp = date_string[:15] + '00' + date_string[17:]
                        date_string = tmp
                    if date_string[18:25] != '000000':
                        seconds = float(date_string[15:17] + '.' + \
                                  date_string[18:25])
                        offset1 = 0.0 - seconds
                        offset2 = 30.0 - seconds
                        offset = min(abs(offset1), abs(offset2))
                    epoch = datetime.datetime.strptime(date_string, \
                            '%y %m %d %H %M %S %f')
                    epochs.append(epoch)
        
        # Check that the RINEX file does not start before the IGS products do,
        # 1 June 1994i
        if good and epochs[0] < igs_start:
            logging.info(f + ': no IGS products available')
            good = False

        # Check that the DOY in file name is right
        if good:
            doy = epochs[0].strftime('%j')
            if f[4:7] != doy:
                logging.info(f + ': incorrect DOY in RINEX file name')
                good = False

        # Calculate the sampling interval
        if good:
            for i in range(1, len(epochs)):
                delta_t = epochs[i] - epochs[i-1]
                deltas.append(delta_t.seconds)
            sampling = statistics.median(deltas)
            if sampling != 30:
                logging.info(f + ': incorrect sampling (' +
                             str(sampling) + ')')
                good = False
        
        # Calculate the duration of the RINEX file
        if good:
            duration = 0
            total_excess = 0
            for delta in deltas:
                if delta >= sampling:
                    duration += sampling
                    excess = delta - sampling
                    total_excess += excess
            # 05:59:00 hours in seconds
            if duration < 21540:
                logging.info(f + ': RINEX duration is less than 6 hours ('
                             + str(duration) + ')')
                good = False
            # 48 hours in seconds    
            if duration > 172800:
                logging.info(f + ': RINEX duration is greater than 48 hours ('
                             + str(duration) + ')')
                good = False
        
        # If the file has passed all the checks write to RinexAntLs.txt
        if good:
            os.rename(f, f.upper())
            stns.add(f.upper()[0:4])
            num_rnx += 1
            fout1.write(f.upper() + ' ' + ant + ' ' + height + '\n')
        
            # Write out: file name, start epoch, end epoch, duration, excess
            l = '{:12s} {} {} {:6d} {:<6d}\n'.\
                format(f.upper(), epochs[0], epochs[-1], int(duration), 
                       int(total_excess))
            fout2.write(l)
        else:
            os.remove(f)

    # Remove if not a RINEX or if the RINEX name is not correctly formed
    elif f not in ['RinexAntLs.txt', 'verifySub.log', 'durations.dat']:
        logging.info(f + ': non-RINEX file or bad NGCA RINEX file name')
        os.remove(f)
fout1.close()
fout2.close()

# Delete verifySub.log if empty
if os.stat('verifySub.log').st_size == 0:
    os.remove('verifySub.log')

# Print stats to the screen
print('There are ' + str(num_rnx) + ' RINEX files and ' + str(len(stns)) +
        ' stations')
