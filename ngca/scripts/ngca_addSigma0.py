#!/usr/bin/env python3

"""Update the Vscale in the NGCA baseline files based on the values in
sigma0.dat
"""

import argparse

from ngca_path import NGCA_DIR

parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")

args = parser.parse_args()

# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'
SIGMA0_FILE = BASELINES_DIR / "sigma0.dat"


# Read in the sigma0 values
sigma = {}
with SIGMA0_FILE.open("r") as f:
    for line in f:
        clus, x, y, sigma0 = line.split()
        clus = clus.replace('_all.simult.adj:Rigorous', '')

        if float(sigma0) < 1:
            sigma[clus] = '1.000'
        else:
            sigma[clus] = sigma0

# Loop over the msr files
for file in BASELINES_DIR.glob('*msr.xml'):

    # Find cluster name
    clus = file.name.replace('_msr.xml', '')
    clus = clus.replace('_all', '')
    
    if clus not in sigma:
        print(f"No sigma0 for {clus}")
        continue

    # Find vscale value for cluster and replace in text
    vscale = sigma[clus]

    text = file.read_text()

    updated_text = text.replace(
        "<Vscale>1.000",
        f"<Vscale>{vscale}",
    )

    file.write_text(text)

    # Write to temp file incase issue with writting
    temp_file = file.with_suffix(".tmp")
    temp_file.write_text(updated_text)

    # Replace original file with tmp file
    temp_file.replace(file)
