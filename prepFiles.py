#!/usr/bin/env python3

"""This script prepares all the RINEX files and antenna information files for
automatic processing
"""
#v2 changes - removing provisions for renaming of files.

import re
import os
import datetime

import sys


# Compile the regular expressions
p1 = re.compile(r'^\w{8}\.\d{2}o$', re.I)

# Create the directory for the individual RINEX antenna information files
os.mkdir('rinexantls')

# Rename files that don't have a '0' before the extension. If renaming will
# cause a name conflict then rename totally
# v2 - Do not rename any input files.
# filenum = 0
# newFile = {}
# fout = open('nameChanges.dat', 'w')
# for f in os.listdir('.'):
#     if p1.match(f):
#         if f[7:8] != '0':
#             new_f = f[:7] + '0' + f[8:]
#             if os.path.isfile(new_f):
#                 conflict = True
#                 while conflict:
#                     new_f = '{:04d}{}0.{}'.format(filenum, f[4:7], f[9:])
#                     filenum += 1
#                     if not os.path.isfile(new_f):
#                         conflict = False
#             newFile[f] = new_f
#             fout.write(f"{new_f} {f}\n")
#             os.rename(f, new_f)
# fout.close()
# this needs to be declared:
newFile = {}
# Read in the start time, end time, and duration for each RINEX file
unsorted_rnx = []
unsorted_start = []
start = {}
end = {}
duration = {}
for line in open('durations.dat'):
    cols = line.split()
    try:
        f = newFile[cols[0]]
    except KeyError:
        f = cols[0]
    unsorted_rnx.append(f)
    try:
        start[f] = datetime.datetime.strptime(cols[1] + ' ' + cols[2], 
                                              '%Y-%m-%d %H:%M:%S')
    except ValueError:
        h, m, s = cols[2].split(':')
        if (abs(float(s) - 0) <= abs(float(s) - 30)):
            s = '00'
        else:
            s = '30'
        cols[2] = h + ':' + m + ':' + s
        start[f] = datetime.datetime.strptime(cols[1] + ' ' + cols[2], 
                                              '%Y-%m-%d %H:%M:%S')
    unsorted_start.append(start[f])
    try:
        end[f] = datetime.datetime.strptime(cols[3] + ' ' + cols[4], 
                                            '%Y-%m-%d %H:%M:%S')
    except ValueError:
        h, m, s = cols[4].split(':')
        if (abs(float(s) - 0) <= abs(float(s) - 30)):
            s = '00'
        else:
            s = '30'
        cols[4] = h + ':' + m + ':' + s
        end[f] = datetime.datetime.strptime(cols[3] + ' ' + cols[4], 
                                            '%Y-%m-%d %H:%M:%S')
    duration[f] = int(float(cols[5]))
sorted_start, sorted_rnx = \
        (list(t) for t in zip(*sorted(zip(unsorted_start, unsorted_rnx))))

# Read in the station, antenna, and height for each RINEX file
station = {}
antenna = {}
height = {}
for line in open('RinexAntLs.txt'):
    cols = line.split()
    try:
        f = newFile[cols[0]]
    except KeyError:
        f = cols[0]
    station[f] = f[:4]
    antenna[f] = cols[1]
    height[f] = cols[2]

# Loop until all RINEX files are in a cluster
while sorted_rnx:
    #print(sorted_rnx[0], sorted_rnx[1], sorted_rnx[2])
    rnx0 = sorted_rnx[0]
    cutoff = end[rnx0]
    data = []
    data.append([rnx0, antenna[rnx0], height[rnx0]])
    sorted_rnx.remove(rnx0)
    #print('rem 1', sorted_rnx[0], sorted_rnx[1], sorted_rnx[2])
        
    # Open the cluster RINEX antemnna information file
    rnxantls = rnx0[9:11] + rnx0[4:7] + '_ls'
    if os.path.isfile(rnxantls):
        conflict = True
        filenum = 0
        while conflict:
            filenum += 1
            rnxantls = rnxantls[:5] + str(filenum) + '_ls'
            if not os.path.isfile(rnxantls):
                conflict = False
    fout = open(rnxantls, 'w')

    # Loop over the RINEX files left
    rnx_remove = []
    for rnx in sorted_rnx:
        #print(sorted_rnx[0], sorted_rnx[1], sorted_rnx[2])
        delta = cutoff - start[rnx]
        #print(rnx, delta.total_seconds()) 
        if delta.total_seconds() >= 7200:
            if end[rnx] < cutoff:
                cutoff = end[rnx]
            data.append([rnx, antenna[rnx], height[rnx]])
            rnx_remove.append(rnx)
            #print('rem 2', sorted_rnx[0], sorted_rnx[1], sorted_rnx[2])
    for rnx in rnx_remove:    
        sorted_rnx.remove(rnx)

    # Write out the cluster data
    for line in data:
        fout.write(line[0] + ' ' + line[1] + ' ' + line[2] + '\n')
    fout.close()

# Move antenna information files to rinexantls/
os.system('mv *_ls rinexantls')

# Delete nameChanges.dat if empty
# v2 - do not check / delete
# if os.stat('nameChanges.dat').st_size == 0:
#     os.remove('nameChanges.dat')
