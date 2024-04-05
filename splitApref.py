#!/usr/bin/env python3

""" This script splits the APREF solution in two. It creates two SINEX files,
one for stations with more than 2 years of data and one for the rest. The second
one keeps two stations from the first as constraints
 
Requirements: 
    - GMT program gmtselect
	- GA program rdsinex

Usage: splitAPREF.py
"""

from geodepy import gnss
import os 


# Read in the solution epoch
epochs = {}
epochs_block = gnss.read_sinex_solution_epochs_block('SNXEPO.SNX')
for epoch in epochs_block:
    if epoch[1:9] != 'SOLUTION' and epoch[0] != '*':
        stat = epoch[1:5]
        start = epoch[16:28]
        if start[0:1] == '9':
            start = '19' + start
        else:
            start = '20' + start
        end = epoch[29:41]
        if end[0:1] == '9':
            end = '19' + end
        else:
            end = '20' + end
        try:
            epochs[stat].append(start)
            epochs[stat].append(end)
        except KeyError:
            epochs[stat] = []
            epochs[stat].append(start)
            epochs[stat].append(end)

# Calculate the time series durations and add appropriate list: con for 2 years
# or more and noncon for the rest
exclude = ['CEDU', 'MAC1', 'HOBA', 'CA19']
con = []
noncon = []
for key, value in epochs.items():
    value.sort()
    duration = int(value[-1][0:4]) - int(value[0][0:4])
    if value[0][5:8] > value[-1][5:8]:
        duration -= 1
    # These two stations link the con and noncon networks, as they are
    # included in both
    if key == 'KALG' or key == '21NA':
        continue       
    if duration < 2 or key in exclude:
        noncon.append(key)
    else:
        con.append(key)

# Create the two SINEX files with the stns in the supplied list removed
gnss.remove_stns_sinex('SNXEPO.SNX', con)
os.rename('output.snx', 'SNXEPO.SNX.NONCON.AUS')
gnss.remove_stns_sinex('SNXEPO.SNX', noncon)
os.rename('output.snx', 'SNXEPO.SNX.CON.AUS')
