#!/usr/bin/bash

# ---
# DESCRIPTION
# - This script will move the SINEX files from 'sinex' folder to
#   the 'solutions' folder in the archive.
# - It will also change the name to the cluster ID format identified
#   in the log_4_submit.txt file.
# - The log_4_submit file should be copied into the 'reupload' folder.
# - It can be used when jobs have failed to trigger and they
#   have been reuploaded (reupload_auspos_job_from_list.sh) and
#   then redownloaded (download_auspos_jobs_from_list.sh).
# - It is aimed to run from the geodesy@[IP ADRESS] ec2.
# - It will move SINEX from jobs listed in list.txt
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
# RUN LOOP TO MOVE SINEX

echo "--- START MOVE"

while read line; do

        # Get job number
        JOB=`echo -e "${line}" | awk '{print $1}'`
        echo "${JOB}"

        # Get cluster ID
        CLUSTER_ID=`grep ${JOB} log_4_submit.txt | awk '{print $6}'`
        ID=${CLUSTER_ID::-3}

        # Move SINEX and change name
        mv ${DIR}/sinex/${JOB}.SNX ../20*/solutions/${ID}.SNX

done < "${INPUT_LIST}"

# ---
# FINISH

echo "--- MOVE COMPLETE"
