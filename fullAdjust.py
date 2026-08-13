#!/usr/bin/env python3

# This script automates the fully constrained adjustment performed after NGCA scaling

# renameless 20251111 - include E N coordinates via stn-coord-types and retain .xyz files for use with DynaDiff

import sys
import os
import glob
import re
import subprocess
from pathlib import Path


def run(cmd):
    """Run a shell command and return True if successful."""
    result = subprocess.run(cmd, shell=True)
    return result.returncode == 0

# ROOT directory (default to ~/GA if not set)
ROOT = Path(os.path.expandvars(os.environ.get("ROOT", str(Path.home() / "GA"))))

# Ensure DynaML.xsd is available in current directory
xsd_source = ROOT / 'aux-files' / 'DynaML.xsd'
xsd_link = Path('DynaML.xsd')
if not xsd_link.exists() and xsd_source.exists():
    xsd_link.symlink_to(xsd_source)

d = os.getcwd()
dir_name = os.path.basename(d)

# Detect Format B: directory named {JUR}_NGCA_{YYYYMMDD}
format_b_match = re.match(r'^([A-Z]{2,3})_NGCA_(\d{8})$', dir_name)

if format_b_match:
    # Format B: extract from directory name
    jur = format_b_match.group(1)
    ngcaVer = format_b_match.group(2)
else:
    # Format A: extract from parent directory and sibling dirs
    jur = d.split('/')[-2].upper()
    p1 = re.compile(r'\d{8}$')
    os.chdir('../')
    dirs = os.listdir('.')
    dirs.sort()
    for dir in dirs:
        if p1.match(dir):
            ngcaVer = dir
    os.chdir(d)

aprefPath = glob.glob(str(ROOT / 'apref' / 'apref*.snx'))
aprefFile = aprefPath[-1]  # Full path to apref SINEX
aprefBasename = Path(aprefFile).name  # e.g., "apref20240713.snx"
discontsFile = str(ROOT / 'apref' / aprefBasename.replace('apref', 'disconts'))

if not run('dnaimport -n ' + jur + '_RENAME ' + jur + '_NGCA_' +
        ngcaVer + '_stn.xml ' + jur + '_NGCA_' + ngcaVer +
        '_msr.xml -r GDA2020 --flag-unused-stations --remove-ignored-msr '
        '--export-xml'):
    sys.exit('ERROR: dnaimport (RENAME) failed')

run('cp ' + jur + '_RENAMEmsr.xml ' + jur + '_GDA2020_' + ngcaVer + '_msr.xml')

if not run('dnaimport -n ' + jur + '_GDA2020 ' + aprefFile +
        ' ' + jur + '_RENAMEstn.xml ' + jur +
        '_RENAMEmsr.xml -r GDA2020 --discontinuity-file ' + discontsFile):
    sys.exit('ERROR: dnaimport (GDA2020) failed')

if not run('dnasegment ' + jur + '_GDA2020 --min 800 --max 800'):
    sys.exit('ERROR: dnasegment failed')

if not run('dnaadjust ' + jur + '_GDA2020 --phased --max-iter 30 '
        '--stn-coord-types "ENzPLHhXYZ" '
        '--output-adj-msr --sort-adj-msr-field 7 --export-xml-stn '
        '--export-xml-msr'):
    sys.exit('ERROR: dnaadjust failed')

if not run('dnaimport -n ' + jur + '_GDA2020_' + ngcaVer +
        ' ' + jur + '_GDA2020.phased.adj.stn.xml ' + jur + '_GDA2020_' +
        ngcaVer + '_msr.xml -r GDA2020 --flag-unused-stations --export-xml'):
    sys.exit('ERROR: dnaimport (final) failed')

run('rm ' + jur + '_RENAME*')
run('cp ' + jur + '_GDA2020.phased.adj temp')
run('cp ' + jur + '_GDA2020.phased.xyz temp2')
run('rm ' + jur + '_GDA2020[.-]*')
run('mv temp ' + jur + '_GDA2020.phased.adj')
run('mv temp2 ' + jur + '_GDA2020.phased.xyz')
run('rm ' + jur + '_GDA2020_' + ngcaVer + '.*')
run('mv ' + jur + '_GDA2020_' + ngcaVer + 'stn.xml ' + jur + '_GDA2020_' +
        ngcaVer + '.adj.xml')
run('rm ' + jur + '_GDA2020_' + ngcaVer + 'msr.xml')
