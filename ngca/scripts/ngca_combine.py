#!/usr/bin/env python3

import argparse
import subprocess

from ngca_path import NGCA_DIR, DYNADJUST_ENV


parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")
parser.add_argument("-a", "--archive", required=True, help="Archive folder to be processed")

args = parser.parse_args()


# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'
ALL_DIR = BASELINES_DIR / "all"

CHECK_SCRIPT = "ngca_check.py"


date = args.archive.strip()

network = f"{args.juris.upper()}_NGCA_{date}"

# Run ngca_check.py
log_file = BASELINES_DIR / "checkNgca.log"

with log_file.open("w") as log:
    subprocess.run(
        [
            str(CHECK_SCRIPT),
            "20",
            "20",
            "40",
            "--directory",
            str(BASELINES_DIR),
        ],
        stdout=log,
        check=True,
    )


# Create the all directory
ALL_DIR.mkdir(exist_ok=True)

# Move *all*.xml files
for file in BASELINES_DIR.glob("*all*.xml"):
    file.rename(ALL_DIR / file.name)

combine_log = BASELINES_DIR / "combine.log"

with combine_log.open("w") as log:

    # Run dnaimport
    subprocess.run(
        [
            "dnaimport",
            "-n",
            network,
            *[str(f) for f in BASELINES_DIR.glob("*stn.xml")],
            *[str(f) for f in BASELINES_DIR.glob("*msr.xml")],
            "-r",
            "ITRF2014",
            "--export-xml",
        ],
        check=True,
        env=DYNADJUST_ENV,
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

    # Run dnareftran
    subprocess.run(
        [
            "dnareftran",
            network,
            "-r",
            "GDA2020",
            "--export-xml",
        ],
        check=True,
        env=DYNADJUST_ENV,
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
    )

# Rename EoO files
(BASELINES_DIR / f"{network}stn.xml").rename(
    BASELINES_DIR / f"{network}_EoO_stn.xml"
)

(BASELINES_DIR / f"{network}msr.xml").rename(
    BASELINES_DIR / f"{network}_EoO_msr.xml"
)

# Rename GDA2020 files
try:
    (BASELINES_DIR / f"{network}.GDA2020stn.xml").rename(
        BASELINES_DIR / f"{network}_stn.xml"
    )

    (BASELINES_DIR / f"{network}.GDA2020msr.xml").rename(
        BASELINES_DIR / f"{network}_msr.xml"
    )
except:
    print("Couldn't rename GDA2020 files")

# Remove remaining network.* files
for file in BASELINES_DIR.glob(f"{network}.*"):
    file.unlink()