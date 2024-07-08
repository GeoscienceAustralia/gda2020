#!/usr/bin/env python3

"""Perform QA on an APREF solution its discontinuity file.

Usage: apref_qa [<date1> <date2>]

Where:
    * date1 is the date of the new APREF solution in YYYYMMDD
    * date2 is the date of the old APREF solution in YYYYMMDD

If two dates are not provided, the most-recent solution will be
checked and the compared with the second-most-recent solution. The
code requires that the necessary files are kept in sub-directories
that are named according to the date of the solution (YYYYMMDD).

This script does some basic checks on the APREF time series
combination solution and associated discontinuity file given by
<date1>. It first checks the files produced by the Positioning Section,
XVSOLFIN.SNX and soln.snx, and then those that are used in the national
adjustment, aprefYYYYMMDD.snx and discontsYYYYMMDD.snx.

The checks include:
    * No duplicate stations, excluding different point codes
    * Matching discontinuities in both files
    * No duplicate DOMES numbers
    * Properly formed epochs
    * Properly formed DOMES numbers

The script also compares the new solution and its discontinuity file (<date1>)
to the old ones (<date2>), reporting on:
    * Stations that have been added or removed from the new solution
    * Stations that have a changed start epoch
    * Stations that have changed DOMES number
    * Discontinuities that have been changed
    * Stations that have moved by more than 1cm in any one direction
    * Positional uncertainties that have changed by +/-10%

The positions of the 109 RVS stations in the files that are used in the
national adjustment are compared with their gazetted values. Stations
that have moved by more than their combined uncertainties (95% CL) are
listed.
"""

import sys
import glob
import os
from numpy import ndarray

from geodepy.gnss import read_solution_epochs
from geodepy.gnss import read_sinex_estimate
from geodepy.gnss import read_sinex_sites
from geodepy.gnss import read_sinex_matrix
from geodepy.gnss import read_disconts
from geodepy.convert import xyz2llh
from geodepy.geodesy import xyz2enu
from geodepy.statistics import vcv_cart2local

# Create the data structures needed to eliminate known duplicate DOMES numbers
duplicate_domes = ['24901S002', '49805S002', '50254S001', '42323S001',
                   '13235S002', '23904S002', '30310S002', '50103M108',
                   '40451S010', '66010M001', '50134M001', '50109S001',
                   '42202M005']
duplicate_domes_stations = {'24901S002': ['BHR3', 'BHR4'],
                            '49805S002': ['EIL3', 'EIL4'],
                            '50254S001': ['MRL1', 'MRL2'],
                            '42323S001': ['MTV1', 'MTV2'],
                            '13235S002': ['OAK1', 'OAK2'],
                            '23904S002': ['OSN3', 'OSN4'],
                            '30310S002': ['PRE3', 'PRE4'],
                            '50103M108': ['TID1', 'TIDB'],
                            '40451S010': ['WDC5', 'WDC6'],
                            '66010M001': ['DAV1', 'DAVR'],
                            '50134M001': ['DARW', 'DARR'],
                            '50109S001': ['ADE1', 'ADE2'],
                            '42202M005': ['AREQ', 'AREV']}

# Determine the dates of the old and the new solution
if len(sys.argv) == 3:
    new = sys.argv[1]
    old = sys.argv[2]
    if new < old:
        print('The date of the old solution is more recent than the date '
              + 'of the new solution')
        sys.exit()
else:
    dirs = glob.glob('20??????')
    dirs.sort()
    new = dirs[-1]
    old = dirs[-2]

# Get the paths to the files
old_sol = None
for f in glob.glob(old + '/XVSOLFIN_*.SNX'):
    old_sol = f
if old_sol is None:
    raise TypeError
old_dis = None
for f in glob.glob(old + '/soln_*.snx'):
    old_dis = f
if old_dis is None:
    raise TypeError
new_sol = None
for f in glob.glob(new + '/AUS0OPSSNX_*_SOL.SNX'):
    new_sol = f
if new_sol is None:
    raise TypeError
new_dis = None
for f in glob.glob(new + '/AUS0OPSSNX_*_DSC.SNX'):
    new_dis = f
if new_dis is None:
    raise TypeError

# Open the output file and write header info
fout = open('xvsolfin' + new + '.txt', 'w')
fout.write('1.0 QA report for XVSOLFIN(' + new + ')\n\n')

# Read in the files
sites_new = read_sinex_sites(new_sol)
epochs_new = read_solution_epochs(new_sol)
estimates_new = read_sinex_estimate(new_sol)
# matrices_new = read_sinex_matrix(new_sol)
disconts_new = read_disconts(new_dis)
sites_old = read_sinex_sites(old_sol)
epochs_old = read_solution_epochs(old_sol)
estimates_old = read_sinex_estimate(old_sol)
# matrices_old = read_sinex_matrix(old_sol)
disconts_old = read_disconts(old_dis)

# Check for duplicate solutions, that is, that each
# 'site code + point code' is unique
num = {}
point_snx = {}
for estimate in estimates_new:
    site = estimate[0]
    point = estimate[1]
    try:
        point_snx[site].append(point)
    except KeyError:
        point_snx[site] = []
        point_snx[site].append(point)
    try:
        num[site + point] += 1
    except KeyError:
        num[site + point] = 1
first = True
for key in num.keys():
    if num[key] > 1:
        if first:
            fout.write(
                '1.1 The following solutions are duplicated:\n')
            first = False
        fout.write('Site: ' + key[0:4] + ' Point: ' + key[4:] + '\n')
if first:
    fout.write('1.1 All solutions are unique\n')
fout.write('\n')

# Check that every discontinuity in the SINEX file is also in the
# discontinuities file
point_dis = []
first = True
for discont in disconts_new:
    site = discont[0]
    point = discont[2]
    point_dis.append(site + point)
for site in point_snx.keys():
    if len(point_snx[site]) == 1 and point_snx[site][0] == '1':
        pass
    else:
        for point in point_snx[site]:
            sol = site + point
            if sol not in point_dis:
                if first:
                    fout.write('1.2 The following solutions are missing from ')
                    fout.write('the discontinuity file:\n')
                    first = False
                fout.write('Site: ' + site + ' Point: ' + point + '\n')
if first:
    fout.write('1.2 All discontinuities are accounted for\n')
fout.write('\n')

# Check for duplicate DOMES numbers
stn = {}
for site in sites_new:
    try:
        stn[site[2]].append(site[0])
    except KeyError:
        stn[site[2]] = [site[0]]
first = True
for key in stn.keys():
    if len(stn[key]) > 1:
        if key in duplicate_domes and \
                stn[key][0] in duplicate_domes_stations[key] and \
                stn[key][1] in duplicate_domes_stations[key]:
            pass
        else:
            if first:
                fout.write('1.3 The following DOMES numbers are used by ')
                fout.write('multiple stations:\n')
                first = False
            dup_stn = ' '.join(stn[key])
            fout.write(key + ': ' + dup_stn + '\n')
if first:
    fout.write('1.3 All DOMES numbers are unique\n')
fout.write('\n')

# Check the format of the DOMES numbers
first = True
for site in sites_new:
    strip_site = site[2].strip()
    if len(strip_site) != 9:
        if first:
            fout.write('1.4 The following DOMES numbers are malformed:\n')
            first = False
        fout.write('Site: ' + site[0] + ' DOMES: ' + site[2])
if first:
    fout.write('1.4 All DOMES numbers are correctly formed')
fout.write('\n\n')

# Check for properly formed solution epochs
first = True
for epoch in epochs_new:
    site = epoch[0]
    sol = epoch[2]
    start = epoch[4]
    end = epoch[5]
    mean = epoch[6]
    if (start == '  :   :0    ' or end == '  :   :0    '
            or end == '  :   :0'):
        if first:
            fout.write('1.5 The following epochs are malformed:\n')
            first = False
        fout.write('Site: ' + site + ' Solution: ' + sol + ' Epochs: "' + start
                   + '", "' + end + '", "' + mean + '"\n')
if first:
    fout.write('1.5 All epochs are correctly formed\n')
fout.write('\n===\n')

fout.write('2.0 Comparison of XVSOLFIN(' + new + ') with XVSOLFIN(' + old)
fout.write(')\n\n')

# List stations that have been added or removed from the new solution
sols_new = set()
sols_old = set()
for estimate in estimates_new:
    sols_new.add(estimate[0] + estimate[1])
for estimate in estimates_old:
    sols_old.add(estimate[0] + estimate[1])
added_sols = list(sols_new - sols_old)
removed_sols = list(sols_old - sols_new)
if added_sols:
    fout.write('2.1 The following solutions have been added:\n')
    for added_sol in added_sols:
        fout.write('Site: ' + added_sol[0:4] + ' Solution: ' + added_sol[4:]
                   + '\n')
else:
    fout.write('2.1 No solutions have been added\n')
fout.write('\n')
if removed_sols:
    fout.write('2.2 The following solutions have been removed:\n')
    for removed_sol in removed_sols:
        fout.write('Site: ' + removed_sol[0:4] + ' Solution: '
                   + removed_sol[4:] + '\n')
else:
    fout.write('2.2 No solutions have been removed\n')
fout.write('\n')

# Check to see if the start epoch of any solutions has changed
data_old = {}
for epoch in epochs_old:
    data_old[epoch[0] + epoch[2]] = {}
    data_old[epoch[0] + epoch[2]]['start'] = epoch[4]
epoch_mismatches = []
for epoch in epochs_new:
    try:
        if epoch[4] != data_old[epoch[0] + epoch[2]]['start']:
            epoch_mismatches.append([epoch[0], epoch[2],
                                     data_old[epoch[0] + epoch[2]]['start'],
                                     epoch[4]])
    except KeyError:
        pass
if epoch_mismatches:
    fout.write('2.3 The start epoch of the following solutions have ')
    fout.write('changed:\n')
    for epoch_mismatch in epoch_mismatches:
        fout.write('Site: ' + epoch_mismatch[0] + ' Solution: '
                   + epoch_mismatch[1] + ' Old: ' + epoch_mismatch[2]
                   + ' New: ' + epoch_mismatch[3] + '\n')
else:
    fout.write('2.3 No change to the start epochs of the solutions\n')
fout.write('\n')

# Check to see if the DOMES number of any site has changed
domes_old = {}
for site in sites_old:
    domes_old[site[0]] = site[2]
dome_mismatches = []
for site in sites_new:
    try:
        if site[2] != domes_old[site[0]]:
            dome_mismatches.append([site[0], domes_old[site[0]], site[2]])
    except KeyError:
        pass
if dome_mismatches:
    fout.write('2.4 The DOMES number of the following sites have changed:\n')
    for dome_mismatch in dome_mismatches:
        fout.write('Site: ' + dome_mismatch[0] + ' Old: ' + dome_mismatch[1]
                   + ' New: ' + dome_mismatch[2] + '\n')
else:
    fout.write('2.4 No changes to any site DOMES number\n')
fout.write('\n')

# Check for discontinuities that have been changed
dis_old = {}
for discont in disconts_old:
    if discont[6] == 'P':
        dis_old[discont[0] + discont[2]] = discont[4] + ' ' + discont[5]
first = True
discont_mismatches = []
for discont in disconts_new:
    if discont[6] == 'P':
        try:
            if discont[4] + ' ' + discont[5] \
                    != dis_old[discont[0] + discont[2]]:
                discont_mismatches.append([discont[0], discont[2],
                                           dis_old[discont[0] + discont[2]],
                                           discont[4], discont[5]])
        except KeyError:
            pass
if discont_mismatches:
    fout.write('2.5 The following discontinuities have been changed:\n')
    for discont_mismatch in discont_mismatches:
        fout.write('Site: ' + discont_mismatch[0] + ' Solution: '
                   + discont_mismatch[1] + ' Old: ' + discont_mismatch[2]
                   + ' New: ' + discont_mismatch[3] + ' '
                   + discont_mismatch[4] + '\n')
else:
    fout.write('2.5 No discontinuities were changed\n')
fout.write('\n')

# Check for stations that have moved by more than 1cm in any one
# direction (enu)
data_old = {}
for estimate in estimates_old:
    data_old[estimate[0] + estimate[1]] = [estimate[3], estimate[4],
                                           estimate[5]]
first = True
for estimate in estimates_new:
    try:
        del_x = estimate[3] - data_old[estimate[0] + estimate[1]][0]
        del_y = estimate[4] - data_old[estimate[0] + estimate[1]][1]
        del_z = estimate[5] - data_old[estimate[0] + estimate[1]][2]
        (lat, lon, h) = xyz2llh(estimate[3], estimate[4], estimate[5])
        (e, n, u) = xyz2enu(lat, lon, del_x, del_y, del_z)
        if abs(e) > 0.1 or abs(n) > 0.1 or abs(u) > 0.1:
            if first:
                fout.write('2.6 The following solutions moved by more than ')
                fout.write('1cm:\n')
                first = False
            line = 'Site: {} Solution: {} Shift (enu): {:.1f} {:.1f} {:.1f}\n' \
                .format(estimate[0], estimate[1], e * 100, n * 100, u * 100)
            fout.write(line)
    except KeyError:
        pass
if first:
    fout.write('2.6 No solutions moved by more than 1cm\n')
fout.write('\n')
"""
# Check for positional uncertainties that have changed by +/-10%
vcv_old = {}
for matrix in matrices_old:
    vcv = ndarray([[float(matrix[2]), float(matrix[3]), float(matrix[4])],
                   float([matrix[3]), float(matrix[5]), float(matrix[6])],
                   float([matrix[4]), float(matrix[6]), float(matrix[7])]])
    vcv_old[matrix[0] + matrix[1]] = vcv
data_old = {}
for estimate in estimates_old:
    (lat, lon, h) = xyz2llh(estimate[3], estimate[4], estimate[5])
    vcv_local = vcv_cart2local(vcv_old[estimate[0] + estimate[1]], lat, lon)
    print(vcv_local)
    data_old[estimate[0] + estimate[1]].append(vcv_local)
    
    
    data_old[estimate[0] + estimate[1]] = [estimate[6], estimate[7],
                                           estimate[8]]
first = True
for estimate in estimates_new:
    try:

        del_x = estimate[3] - data_old[estimate[0] + estimate[1]][0]
        del_y = estimate[4] - data_old[estimate[0] + estimate[1]][1]
        del_z = estimate[5] - data_old[estimate[0] + estimate[1]][2]
        (lat, lon, h) = xyz2llh(estimate[3], estimate[4], estimate[5])
        (e, n, u) = xyz2enu(lat, lon, del_x, del_y, del_z)
        if abs(e) > 0.1 or abs(n) > 0.1 or abs(u) > 0.1:
            if first:
                fout.write('The following uncertainties have changed ')
                fout.write('by +/-10%:\n')
                first = False
            line = 'Site: {} Solution: {} Shift (enu): {:.1f} {:.1f} {:.1f}\n'\
                .format(estimate[0], estimate[1], e*100, n*100, u*100)
            fout.write(line)
    except KeyError:
        pass
"""
"""
    * RVS stations that have moved significantly from their gazettal
        values
"""
