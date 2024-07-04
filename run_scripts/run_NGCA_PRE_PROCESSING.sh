#!/usr/bin/bash

# Archive
ARCHIVE=20240515

# List of jurisdictions to loop through (largest to smallest)
JURIS_LIST=("act" "tas" "vic" "sa" "nt" "qld" "nsw")
#JURIS_LIST=("wa")

# Create automatic notes file (this will remove existing file)
rm -f /home/fedora/ngca/auto_notes_$ARCHIVE.txt
touch /home/fedora/ngca/auto_notes_$ARCHIVE.txt

# Commentary 
printf "\nRunning script to organise NGCA files up until the stage before AUSPOS processing."
printf "\nWill be working on these jurisdictions: ${JURIS_LIST[*]^^}, for the $ARCHIVE NGCA archive."

for JURIS in ${JURIS_LIST[*]}; do

    # Navigate to jurisdiction directory
    printf "\nNavigating to ${JURIS^^} archive directory...\n"
    cd /home/fedora/ngca/$JURIS/$ARCHIVE

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
    echo "${JURIS^^} has $NumRINEX RINEX and $NumJobs jobs" >> /home/fedora/ngca/auto_notes_$ARCHIVE.txt

    # Commentary
    printf "\nFinished working on ${JURIS^^}\n"

done

printf "\n------ DONE ------\n"
printf "\n"

