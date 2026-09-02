#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -d "$SCRIPT_DIR/scripts" ]] && [[ ":$PATH:" != *":$SCRIPT_DIR/scripts:"* ]]; then
    export PATH="$PATH:$SCRIPT_DIR/scripts"
fi

# List of jurisdictions to run consecutively - increasing number of RINEX
#JURIS_LIST=("tas" "act" "vic" "sa" "nt" "wa" "nsw" "qld")

# Parse jurisdiction flag
JURIS_LIST=()

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--juris)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                JURIS_LIST+=("${1,,}")   # force lowercase
                shift
            done
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 -j tas act nsw"
            exit 1
            ;;
    esac
done

# Use all jurisdictions if none specified
if [[ ${#JURIS_LIST[@]} -eq 0 ]]; then
    JURIS_LIST=("tas" "act" "vic" "sa" "nt" "wa" "nsw" "qld")
fi

# Commentary 
printf "\nRunning script to download and organise NGCA files before AUSPOS processing."
printf "\nWill be working on these jurisdictions: ${JURIS_LIST[*]^^}, for the $ARCHIVE NGCA archive."

printf "\nDownloading NGCA files into the NGCA archive.\n"

while IFS= read -r line; do
    if [[ -z "${ARCHIVE:-}" ]]; then
        ARCHIVE="$line"
    else
        printf '%s\n' "$line"
    fi
done < <(PYTHONUNBUFFERED=1 ngca_get.py -j "${JURIS_LIST[@]}")

# Create automatic notes file (removes existing note file - for a full auto_note summary run all juris consecutively)
rm -f $SCRIPT_DIR/auto_notes_$ARCHIVE.txt
touch $SCRIPT_DIR/auto_notes_$ARCHIVE.txt

for JURIS in ${JURIS_LIST[*]}; do

    # Navigate to jurisdiction directory
    printf "\nNavigating to ${JURIS^^} archive directory...\n"
    cd $SCRIPT_DIR/$JURIS/$ARCHIVE

    # Run verifySub.py
    printf "\nexecuting verifySub.py for ${JURIS^^}...\n"
    ngca_verifySub.py

    # Run prepFiles.py
    printf "\nexecuting prepFiles.py for ${JURIS^^}...\n"
    ngca_prepFiles.py

    # Identify number of jobs
    NumJobs=$(ls rinexantls/ | wc -l)    
    printf "\nNumber of jobs for ${JURIS^^}: $NumJobs.\n"

    # Identify number of RINEX
    NumRINEX=$(ls *.*O | wc -l)

    # Write number of jobs to a notes file"
    echo "${JURIS^^} has $NumRINEX RINEX and $NumJobs jobs" >> $SCRIPT_DIR/auto_notes_$ARCHIVE.txt

    # Commentary
    printf "\nFinished working on ${JURIS^^}\n"

done

printf "\n------ DONE ------\n"
printf "\n"

