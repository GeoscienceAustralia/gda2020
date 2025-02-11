#!/usr/bin/env python3

"""Compare the latest RVS station coordinates to those that were gazetted
Usage: checkRVS.py [<date>]
Where <date> is the date of the APREF solution to be checked in YYYYMMDD
The comparison is done in Cartesian space but the coordinate differences are
given in local space
If a date is not provided, the most-recent solution will be used. The
code requires that the necessary files are kept in sub-directories
that are named according to the date of the solution (YYYYMMDD) and that the
original RVS SINEX file is in a sub-directory named 'rvs'.
"""

import sys
import glob
from geodepy.gnss import read_sinex_estimate
from geodepy.convert import xyz2llh
from geodepy.geodesy import xyz2enu
from math import sqrt


# The date of the APREF cumulative solution to be checked can be passed at the
# command line. If not, determine what is the most recent sub-directory
if len(sys.argv) == 2:
    sol_date = sys.argv[1]
else:
    dirs = glob.glob('20??????')
    dirs.sort()
    sol_date = dirs[-1]

# Create paths to the two input SINEX files
rvs_sol = 'rvs/apref16282_typeb.SNX'
apref_sol = sol_date + '/apref' + sol_date + '.snx'

# Read in the coordinate estimates and standard deviations from both the RVS
# and APREF solutions
rvs_estimates = read_sinex_estimate(rvs_sol)
apref_estimates = read_sinex_estimate(apref_sol)

# Open the output file
fout = open('apref' + sol_date + '.rvs', 'w')
fout.write('Comparison of APREF(' + apref_sol[0:8] +
           ') to the RVS gazettal\n\n')

# Determine the 109 RVS stations (4-character ID + solution number)
rvs_stns = set()
rvs_solns = []
for line in open('rvs/RVS_GDA2020.txt'):
    if line[0] != '#':
        rvs_stns.add(line[0:4])
        rvs_solns.append(line[0:4] + line[6:7])

# Loop over the RVS estimates and add the coordinates and uncertainties to a
# dictionary. The standard deviations are expanded to 95% using the factor of
# 2 suggested by the NMI
rvs_data = {}
for estimate in rvs_estimates:
    if estimate[0] + estimate[1] in rvs_solns:
        rvs_data[estimate[0]] = {'X': estimate[3], 'Y': estimate[4],
                                 'Z': estimate[5], 'sigma_X': 2 * estimate[6],
                                 'sigma_Y': 2 * estimate[7],
                                 'sigma_Z': 2 * estimate[8]}

# Get the number of the last solution for each RVS station in the APREF
# cumulative solution
solns = {}
for estimate in apref_estimates:
    if estimate[0] in rvs_data.keys():
        try:
            solns[estimate[0]].append(estimate[1])
        except KeyError:
            solns[estimate[0]] = []
            solns[estimate[0]].append(estimate[1])
last_soln = []
for key in solns.keys():
    last_soln.append(key + max(solns[key]))

# Loop over the APREF estimates and add the coordinates and uncertainties to a
# dictionary. The standard deviations are expanded to 95% using the factor of
# 2 suggested by the NMI
apref_data = {}
apref_stns = set()
for estimate in apref_estimates:
    if estimate[0] + estimate[1] in last_soln:
        apref_stns.add(estimate[0])
        apref_data[estimate[0]] = {'X': estimate[3], 'Y': estimate[4],
                                   'Z': estimate[5],
                                   'sigma_X': 2 * estimate[6],
                                   'sigma_Y': 2 * estimate[7],
                                   'sigma_Z': 2 * estimate[8]}

# Write out the RVS stations that are missing from the APREF solution
missing = rvs_stns - apref_stns
missing = ', '.join(missing)
if missing:
    fout.write('Missing stations: ' + missing + '\n\n')

# Write out the headers
line = '{} {:>6} {:>5} {:>6} {:>6} {:>6} {:>7} {:>5} {} {:>5} {} {:>5} {}\n'\
    .format('Stn', 'D_3D', 'D_2D', 'East', 'North', 'Up', 'Out?', 'delX',
            'sigX', 'delY', 'sigY', 'delZ', 'sigZ')
fout.write(line)

# Loop over the APREF stations and create a list of dictionaries containing
# the info on each station
data = []
for key in apref_data.keys():

    # Calculate the the coordinate differences in mm
    delX = rvs_data[key]['X'] - apref_data[key]['X']
    delY = rvs_data[key]['Y'] - apref_data[key]['Y']
    delZ = rvs_data[key]['Z'] - apref_data[key]['Z']

    # Calculate the combined uncertainty
    sigmaX = sqrt(rvs_data[key]['sigma_X']**2 + apref_data[key]['sigma_X']**2)
    sigmaY = sqrt(rvs_data[key]['sigma_Y']**2 + apref_data[key]['sigma_Y']**2)
    sigmaZ = sqrt(rvs_data[key]['sigma_Z']**2 + apref_data[key]['sigma_Z']**2)

    # Determine whether the difference is significant or not
    if abs(delX) > abs(sigmaX):
        x_out = 'Y'
    else:
        x_out = 'N'
    if abs(delY) > abs(sigmaY):
        y_out = 'Y'
    else:
        y_out = 'N'
    if abs(delZ) > abs(sigmaZ):
        z_out = 'Y'
    else:
        z_out = 'N'

    # Calculate lat/lon from XYZ
    geographic = xyz2llh(rvs_data[key]['X'], rvs_data[key]['Y'],
                         rvs_data[key]['Z'])
    lat = geographic[0]
    lon = geographic[1]

    # Convert from XYZ to enu (in mm)
    local = xyz2enu(lat, lon, delX, delY, delZ)
    east = local[0] * 1000
    north = local[1] * 1000
    up = local[2] * 1000

    # Calculate the distance and the on ground distance
    dist = sqrt(east**2 + north**2 + up**2)
    on_ground = sqrt(east**2 + north**2)

    # Create a dictionary for this station (converting everything to mm)
    info = {'stn': key, 'dist': dist, 'on_ground': on_ground, 'east': east,
            'north': north, 'up': up, 'x_out': x_out, 'y_out': y_out,
            'z_out': z_out, 'delX': delX, 'sigmaX': sigmaX, 'delY': delY,
            'sigmaY': sigmaY, 'delZ': delZ, 'sigmaZ': sigmaZ}
    data.append(info)

# Sort the list of dictionaries, format and write to output file
data.sort(key=lambda x: x.get('dist'), reverse=True)
for info in data:
    line = '{} {:5.1f} {:5.1f} {:6.1f} {:6.1f} {:6.1f} {:>3} {} {} {:5.1f} ' \
           '{:4.1f} {:5.1f} {:4.1f} {:5.1f} {:4.1f}\n'.\
        format(info['stn'], info['dist'], info['on_ground'], info['east'],
               info['north'], info['up'], info['x_out'], info['y_out'],
               info['z_out'], info['delX']*1000, info['sigmaX']*1000,
               info['delY']*1000, info['sigmaY']*1000, info['delZ']*1000,
               info['sigmaZ']*1000)
    fout.write(line)