#!/usr/bin/env python3

"""Package NGCA results and upload them to the AGRS Portal."""

import argparse
import shutil
import subprocess
from pathlib import Path

from ngca_path import NGCA_DIR, NGCA_SENT_DIR, NGCA_TRANS_TABLES_DIR


parser = argparse.ArgumentParser(
    description="Package NGCA results and upload them to the AGRS Portal.",
)

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")
parser.add_argument("-a", "--archive", required=True, help="Archive folder to be processed")

args = parser.parse_args()

jur = args.juris.lower()
jur_upper = jur.upper()
date = args.archive.strip()

# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'
SINEX_DIR = NGCA_DIR / args.juris.lower() / 'sinexFiles'
RINEXANTLS_DIR = NGCA_DIR / args.juris.lower() / 'rinexantls'


if not NGCA_SENT_DIR.is_dir():
    raise FileNotFoundError(
        f"Sent directory does not exist: {NGCA_SENT_DIR}"
    )


DATE_DIR = JURIS_DIR / date

campaign = f"{jur_upper}_NGCA_{date}"

CAMPAIGN_DIR = NGCA_SENT_DIR / campaign
ZIP_FILE = NGCA_SENT_DIR / f"{campaign}.zip"


# Campaign subdirectories
DYNAML_DIR = CAMPAIGN_DIR / "dynaML"
ADJUSTMENTS_DIR = CAMPAIGN_DIR / "adjustments"
FEEDBACK_DIR = CAMPAIGN_DIR / "feedback"
OTHER_DIR = CAMPAIGN_DIR / "other"
XYZ_DIR = CAMPAIGN_DIR / "xyz"
CAMPAIGN_SNX_DIR = CAMPAIGN_DIR / "snx"
CLUSTERS_DIR = CAMPAIGN_DIR / "clusters"

# S3 file location
#S3_DESTINATION = "s3://gda2020-ngca/2025_Testing/Portal/"
S3_DESTINATION = "s3://gda2020-ngca/2026_Testing/NGCA/"
dry_run = True

print(f"Preparing NGCA: {campaign}")

# Create campaign directory
if CAMPAIGN_DIR.exists():
    print(f"Removing existing campaign directory: {CAMPAIGN_DIR}")
    shutil.rmtree(CAMPAIGN_DIR)

CAMPAIGN_DIR.mkdir(parents=True)

for directory in (
    DYNAML_DIR,
    ADJUSTMENTS_DIR,
    FEEDBACK_DIR,
    OTHER_DIR,
    XYZ_DIR,
):
    directory.mkdir()


# Copy over xml files

# Copy all XML files into dynaML first
for source in BASELINES_DIR.glob("*.xml"):
    shutil.copy2(source, DYNAML_DIR)


# Remove GDA2020 XML files from dynaML
for file in DYNAML_DIR.glob("*GDA2020*.xml"):
    file.unlink()


# Move NGCA XML files from dynaML into the campaign root
for file in DYNAML_DIR.glob("*NGCA*.xml"):
    destination = CAMPAIGN_DIR / file.name
    file.rename(destination)


# Copy over adjustment files

for file in BASELINES_DIR.glob("*.adj"):
    shutil.copy2(file, ADJUSTMENTS_DIR)


# Copy feedback files

for name in (
    "sigma0.dat",
    "checkNgca.log",
):
    file = BASELINES_DIR / name

    if source.is_file():
        shutil.copy2(file, FEEDBACK_DIR)
    else:
        print(f"Warning: feedback file not found: {source}")


for name in (
    "log_4_submit.txt",
    "selectRef.log",
    "verifySub.log",
):
    source = DATE_DIR / name

    if source.is_file():
        shutil.copy2(source, FEEDBACK_DIR)
    else:
        print(f"Warning: optional feedback file not found: {source}")


# Copy transformation table

trans_tables = sorted(
    NGCA_TRANS_TABLES_DIR.glob(f"{jur}TransTable*.csv")
)

if not trans_tables:
    print(
        "Warning: no transformation tables found matching "
        f"{jur}TransTable*.csv in {NGCA_TRANS_TABLES_DIR}"
    )

for file in trans_tables:
    shutil.copy2(file, OTHER_DIR)


# Copy over xyz files

xyz_files = sorted(BASELINES_DIR.glob("*.xyz"))

if not xyz_files:
    print(f"Warning: no XYZ files found in {BASELINES_DIR}")

for file in xyz_files:
    shutil.copy2(file, XYZ_DIR)


# Copy over sinex files

shutil.copytree(
    SINEX_DIR,
    CAMPAIGN_SNX_DIR,
)


# Match the original "rm snx/*.SNX"
for file in CAMPAIGN_SNX_DIR.glob("*.SNX"):
    file.unlink()
for file in CAMPAIGN_SNX_DIR.glob("*.SNX.ITRF2014"):
    file.unlink()

# Copy over cluster files

if not RINEXANTLS_DIR.is_dir():
    raise FileNotFoundError(
        f"RINEX antenna directory does not exist: "
        f"{RINEXANTLS_DIR}"
    )

shutil.copytree(
    RINEXANTLS_DIR,
    CLUSTERS_DIR,
)


print(f"Creating ZIP archive: {ZIP_FILE.name}")


if ZIP_FILE.exists():
    ZIP_FILE.unlink()

created_archive = shutil.make_archive(
    base_name=str(CAMPAIGN_DIR),
    format="zip",
    root_dir=NGCA_SENT_DIR,
    base_dir=campaign,
)

created_archive = Path(created_archive)

if not created_archive.is_file():
    raise RuntimeError(
        f"ZIP archive was not created: {created_archive}"
    )


# Remove the unpacked campaign only after the ZIP has been created
shutil.rmtree(CAMPAIGN_DIR)

print(f"Created: {created_archive}")

# First remove old file
if not dry_run:
    subprocess.run(
        [
            "aws",
            "s3",
            "rm",
            S3_DESTINATION,
            "--quiet",
            "--recursive",
            "--exclude",
            "*",
            "--include",
            f"*{jur_upper}*",
        ],
        check=True,
    )
else:
    print(f"Would have deleted old file")

print(f"Uploading {created_archive.name}...")

subprocess.run(
    [
        "aws",
        "s3",
        "cp",
        str(created_archive),
        S3_DESTINATION,
        "--quiet",
    ],
    check=True,
)
