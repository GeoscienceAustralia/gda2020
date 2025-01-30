#!/usr/bin/env python3

"""Check the jurisdicitonal translation tables for duplicates"""

import glob
import os
import sys
from pathlib import Path

# Create lists of stations and 4-character IDs to ignore
ignoreStns = ['6635/ 1199', '6626/ 1035', '6429/ 2959', '6731/ 1089',
        '6326/ 1586', '5932/ 1315', '6230/ 1234', '6627/22981', '6532/ 1865',
        '6030/ 1200', '5831/ 1565', '6626/ 3925', '6331/ 1393', '5831/ 1291',
        '6628/57113', '6728/ 7304', '6432/ 2690', '6629/ 3826', '6430/ 2314',
        '6628/39005', '6925/ 2549', '5833/ 1244', '6824/ 1439', '6627/24959',
        '6923/ 1784', '6628/59058', '6630/ 3241', '6534/ 1272', '6632/ 1601',
        '6632/ 1600', '6632/ 1602', '6732/ 1286', '6433/ 2322', '6030/ 1625',
        '6531/ 3264', '6727/ 7103', '6532/ 1927', '6526/ 1390', '6630/ 1507',
        '6736/ 1189', '6332/ 1346', '6326/ 1869', '6031/ 1307', '6729/ 1002',
        '6737/ 1117', '6729/ 1169', '5929/ 1273', '6635/ 1198', '6526/ 2025',
        '6627/ 1609', '6630/ 1109', '6631/ 2026', '6631/ 1121', '6635/ 1197',
        '6634/ 1237', '6629/ 1276', '6630/ 1169', '6029/ 1923', '6429/ 1200',
        '6226/ 1417', '6826/ 2413', '6326/ 1619', '6633/ 1145', '6530/ 2110',
        '6627/24960', '6530/ 2115', '6530/ 2114', '6331/ 1394', '6326/ 1870',
        '6327/ 1367', '6729/ 2210', '6729/ 2211', '6727/ 7312', '6130/ 1436',
        '5932/ 1098', '6536/ 1437', '6028/ 3974', '6527/ 8844', '6628/57240',
        '6131/ 1563', '6432/ 1526', '6628/57112', '6432/ 1524', '6838/ 1138',
        '6729/ 1217', '6226/ 1311', '6226/ 1312', '6533/ 1467', '6729/ 2209',
        '6328/ 1594', '6530/ 1202', '6327/ 1247', '6631/ 2024', '6634/ 1236',
        '6631/ 1162', '6631/ 1167', '6426/ 1641', '6526/ 2041', '6637/ 1122',
        '6728/ 1084', '6230/ 1602', '6531/ 1683', '6728/ 1088', '6331/ 1377',
        '6728/ 7303', '6428/ 2095', '5928/ 1370', '6728/ 7104', '6728/ 7101',
        '6728/ 7102', '6432/ 2635', '6531/ 3358', '6730/ 1428', '5929/ 1347',
        '6635/ 1034', '5830/ 1167', '6630/ 3100', '6630/ 3101', '6428/ 2213',
        '6428/ 2214', '6426/ 1319', '6428/ 1097', '6631/ 2025', '6629/ 7911',
        '6433/ 2350', '6728/ 1000']

# Read in the APREF station names
for file in glob.glob(f'{Path.home()}/apref/apref*.disconts'):
    aprefFile = file
aprefStns = set()
for line in open(aprefFile):
    if line[0:1] != '#':
        aprefStns.add(line[0:4])

# Change to the translation tables directory and loop over the translation
# tables
allIds = {}
allStns = {}
trans1 = {}
trans2 = {}
os.chdir(f'{Path.home()}/transTables')
for file in glob.glob('*.csv'):
    print('Checking ' + file + '...')
    state = file[0:2]
    if state == 'vi':
        state = 'vic'
    elif state == 'ns':
        state = 'nsw'
    elif state == 'ac':
        state = 'act'
    elif state == 'ta':
        state = 'tas'
    elif state == 'ql':
        state = 'qld'
    trans1[state] = {}    
    trans2[state] = {}    
    ids = {}
    stns = {}
    for line in open(file):
        line = line.strip()
        if line != '':
            (id, stn) = line.split(',')
            id = id.strip()
            if id in aprefStns:
                if id != stn:
                    if state != 'sa' and id != 'COLL':
                        print(id + ' is an APREF station') 
            stn = stn.strip()
            if stn in aprefStns:
                if id != stn:
                    print(stn + ' is an APREF station') 
            try:
                ids[id].append(stn)
            except KeyError:
                ids[id] = []
                ids[id].append(stn)
            try:
                stns[stn].append(id)
            except KeyError:
                stns[stn] = []
                stns[stn].append(id)
            try:
                allIds[id].append(state)
            except KeyError:
                allIds[id] = []
                allIds[id].append(state)
            try:
                allStns[stn].append(state)
            except KeyError:
                allStns[stn] = []
                allStns[stn].append(state)
            trans1[state][id] = stn
            trans2[state][stn] = id
    for key, value in ids.items():
        if len(value) > 1:
            str = ' '.join(value)
            print(key + ': ' + str)
    for key, value in stns.items():
        if len(value) > 1:
            if key not in ignoreStns:
                str = ' '.join(value)
                print(key + ': ' + str)
print('Checking as a whole...')
good = True
for key, values in allIds.items():
    if len(set(values)) > 1:
        check = set()
        for value in values:
            check.add(trans1[value][key])
        if len(check) > 1:    
            str = ' '.join(values)
            print(key + ': ' + str)
            good = False
for key, values in allStns.items():
    if len(set(values)) > 1:
        check = set()
        for value in values:
            check.add(trans2[value][key])
        if len(check) > 1:
            str = ' '.join(value)
            print(key + ': ' + str)
            good = False
print('All good!')
