#!/usr/bin/env python3

"""Perform the fully constrained adjustment after NGCA scaling."""

import argparse
import subprocess

from ngca_path import NGCA_DIR, NGCA_APREF_DIR, DYNADJUST_ENV


parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")
parser.add_argument("-a", "--archive", required=True, help="Archive folder to be processed")

args = parser.parse_args()


# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'

jur = args.juris.upper()
date = args.archive.strip()

apref_files = sorted(NGCA_APREF_DIR.glob("apref*.snx"))
discont_files = sorted(NGCA_APREF_DIR.glob("disconts*.snx"))

if not apref_files:
    raise FileNotFoundError(f"No apref*.snx found in {JURIS_DIR}")
if not discont_files:
    raise FileNotFoundError(f"No disconts*.snx found in {JURIS_DIR}")

apref_file = apref_files[-1]
discont_file = discont_files[-1]

# Network names
ngca_network = f"{jur}_NGCA_{date}"
rename_network = f"{jur}_RENAME"
gda_network = f"{jur}_GDA2020"
final_network = f"{jur}_GDA2020_{date}"

# Keep DynAdjust output out of the console
log_file = BASELINES_DIR / "fullAdjust.log"

with log_file.open("w") as log:

    subprocess.run(
        [
            "dnaimport",
            "-n",
            rename_network,
            f"{ngca_network}_stn.xml",
            f"{ngca_network}_msr.xml",
            "-r",
            "GDA2020",
            "--flag-unused-stations",
            "--remove-ignored-msr",
            "--export-xml",
        ],
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=DYNADJUST_ENV,
    )


    # Keep copy of measurement file
    source = BASELINES_DIR / f"{rename_network}msr.xml"
    destination = BASELINES_DIR / f"{final_network}_msr.xml"

    destination.write_bytes(source.read_bytes())


    subprocess.run(
        [
            "dnaimport",
            "-n",
            gda_network,
            str(apref_file),
            f"{rename_network}stn.xml",
            f"{rename_network}msr.xml",
            "-r",
            "GDA2020",
            "--discontinuity-file",
            str(discont_file),
        ],
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=DYNADJUST_ENV,
    )


    subprocess.run(
        [
            "dnasegment",
            gda_network,
            "--min",
            "800",
            "--max",
            "800",
        ],
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=DYNADJUST_ENV
    )


    subprocess.run(
        [
            "dnaadjust",
            gda_network,
            "--phased",
            "--max-iter",
            "30",
            "--stn-coord-types",
            "ENzPLHhXYZ",
            "--output-adj-msr",
            "--sort-adj-msr-field",
            "7",
            "--export-xml-stn",
            "--export-xml-msr",
        ],
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=DYNADJUST_ENV
    )


    subprocess.run(
        [
            "dnaimport",
            "-n",
            final_network,
            f"{gda_network}.phased.adj.stn.xml",
            f"{final_network}_msr.xml",
            "-r",
            "GDA2020",
            "--flag-unused-stations",
            "--export-xml",
        ],
        cwd=BASELINES_DIR,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=DYNADJUST_ENV,
    )


# Remove temporary RENAME files
for file in BASELINES_DIR.glob(f"{rename_network}*"):
    file.unlink()

# Keep the final adjustment and XYZ files
adj_file = BASELINES_DIR / f"{gda_network}.phased.adj"
xyz_file = BASELINES_DIR / f"{gda_network}.phased.xyz"
msr_file = BASELINES_DIR / f"{gda_network}.phased.adj.msr.xml"
stn_file = BASELINES_DIR / f"{gda_network}.phased.adj.stn.xml"

adj_data = adj_file.read_bytes()
xyz_data = xyz_file.read_bytes()
msr_data = msr_file.read_bytes()
stn_data = stn_file.read_bytes()

# Remove intermediate GDA2020 files
for pattern in (
    f"{gda_network}.*",
    f"{gda_network}-*",
):
    for file in BASELINES_DIR.glob(pattern):
        file.unlink()

# Restore retained files
adj_file.write_bytes(adj_data)
xyz_file.write_bytes(xyz_data)
msr_file.write_bytes(msr_data)
stn_file.write_bytes(stn_data)

# Remove intermediate final-network files
for file in BASELINES_DIR.glob(f"{final_network}.*"):
    file.unlink()
