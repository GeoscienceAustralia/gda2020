#!/usr/bin/bash

# ---
# DESCRIPTION
# - This script will download the SINEX files from AUSPOS s3.
# - It can be used when jobs have failed to trigger and they
#   have been reuploaded (reupload_auspos_job_from_list.sh).
# - It will download SINEX from jobs listed in list.txt
# - It is aimed to run from the geodesy@[IP ADRESS] ec2.
# - Format is one column for the list.txt file.
#
#   Example:
#
#       0000000476671-00000000
#       0000000476813-00000000

# ---
# SETUP

DIR=/data/craig/reupload
INPUT_LIST=${DIR}/list.txt
cd ${DIR}

# ---
# RUN LOOP TO DOWNLOAD JOBS

echo "--- START DOWNLOAD"

while read line; do

        # Get job number
        JOB=`echo -e "${line}" | awk '{print $1}'`
        echo "${JOB}"

        # Download the file from s3
        aws s3 cp s3://auspos-test-prod/uploads/${JOB}/ ${DIR}/sinex --recursive  --exclude "*" --include "*-*.SNX"

done < "${INPUT_LIST}"

# ---
# FINISH

echo "--- DOWNLOAD COMPLETE"
