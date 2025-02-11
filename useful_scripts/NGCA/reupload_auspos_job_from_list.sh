#!/usr/bin/bash

# ---
# DESCRIPTION
# - This script will reupload the auspos_job.json file to the AUSPOS s3.
# - It is used when jobs have failed to trigger.
# - It will run the jobs listed in list.txt.
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
# RUN LOOP TO REUPLOAD JOBS

echo "--- START REUPLOAD"

while read line; do

        # Get job number
        JOB=`echo -e "${line}" | awk '{print $1}'`
        echo "${JOB}"

        # Download the file from s3
        aws s3 cp s3://auspos-test-prod/uploads/${JOB}/auspos_job.json ${DIR}/auspos_job.json

        # Upload the file back to s3 (to trigger AUSPOS job)
        aws s3 cp ${DIR}/auspos_job.json s3://auspos-test-prod/uploads/${JOB}/auspos_job.json

done < "${INPUT_LIST}"

# ---
# FINISH

rm auspos_job.json

echo "--- REUPLOAD COMPLETE"
