#!/usr/bin/bash

# Archive date for NGCA processing, e.g. the date on which get_ngca.py has downloaded the NGCA archive contents for the jurisdictions being processed 
ARCHIVE=20240715

# List of jurisdictions to run consecutively - increasing number of RINEX
JURIS_LIST=("tas" "act" "vic" "sa" "nt" "wa" "nsw" "qld")

# Alternatively run a single jurisdiction
#JURIS_LIST=("tas")

# Create automatic notes file (removes existing note file - for a full auto_note summary run all juris consecutively)
rm -f ~/ngca/auto_notes_$ARCHIVE.txt
touch ~/ngca/auto_notes_$ARCHIVE.txt

# Commentary 
printf "\nRunning script to organise NGCA files up until the stage before AUSPOS processing."
printf "\nWill be working on these jurisdictions: ${JURIS_LIST[*]^^}, for the $ARCHIVE NGCA archive."

for JURIS in ${JURIS_LIST[*]}; do

    # Navigate to jurisdiction directory
    printf "\nNavigating to ${JURIS^^} archive directory...\n"
    cd ~/ngca/$JURIS/$ARCHIVE

    # Run verifySub.py
    printf "\nexecuting verifySub.py for ${JURIS^^}...\n"
    verifySub.py

    # Run prepFiles.py
    printf "\nexecuting prepFiles.py for ${JURIS^^}...\n"
    prepFiles.py

    # Identify number of jobs
    NumJobs=$(ls rinexantls/ | wc -l)    
    printf "\nNumber of jobs for ${JURIS^^}: $NumJobs.\n"

    # Identify number of RINEX
    NumRINEX=$(ls *.*O | wc -l)

    # Write number of jobs to a notes file"
    echo "${JURIS^^} has $NumRINEX RINEX and $NumJobs jobs" >> ~/ngca/auto_notes_$ARCHIVE.txt

    # Commentary
    printf "\nFinished working on ${JURIS^^}\n"

done

printf "\n------ DONE ------\n"
printf "\n"

