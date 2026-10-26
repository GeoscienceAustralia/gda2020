#!/usr/bin/env python3
  
"""This script updates the translation tables in transTables/ with those in the
S3 bucket
"""

from ngca_path import NGCA_TRANS_TABLES_DIR
import subprocess

# Define the jurisdictions
jurisdictions = ['nsw', 'act', 'vic', 'tas', 'sa', 'wa', 'nt', 'qld']

# Remove all file in the transTables directory
for item in NGCA_TRANS_TABLES_DIR.iterdir():
    if item.is_file():
        item.unlink()
        
# Copy the files directory from the S3 bucket
subprocess.run(['aws', 's3', 'cp', 's3://gda2020-ngca/files/', str(NGCA_TRANS_TABLES_DIR), '--recursive', '--exclude', '*', '--include', '*Trans*'], check=True)

# Loop over the csv files
for item in NGCA_TRANS_TABLES_DIR.iterdir():

    f = item.name

# Determine the jurisdiction
    if f[:3].lower() in jurisdictions:
        jurisdiction = f[:3].lower()
    elif f[:2].lower() in jurisdictions:
        jurisdiction = f[:2].lower()
    else:
        print(f'File {f} has an unknown jurisdiction')
        item.unlink()

    # Determine the date
    date = f.split('.')[0]
    date = date[-8:]

    # Create the new file name 
    f_new = NGCA_TRANS_TABLES_DIR / f'{jurisdiction}TransTable{date}.csv'
    item.rename(f_new)
