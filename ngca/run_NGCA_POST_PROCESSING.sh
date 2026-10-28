#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -d "$SCRIPT_DIR/scripts" ]] && [[ ":$PATH:" != *":$SCRIPT_DIR/scripts:"* ]]; then
    export PATH="$PATH:$SCRIPT_DIR/scripts"
fi

# TEMP
CPUS=4

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
        -a|--archive)
            shift
            ARCHIVE="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 -j tas act nsw -a 20260825"
            exit 1
            ;;
    esac
done

# Check if archive date is specified
if [[ -z "$ARCHIVE" ]]; then
    echo "Archive date not specified. Use -a <archive_date> to specify the archive date."
    exit 1
fi

# Use all jurisdictions if none specified
if [[ ${#JURIS_LIST[@]} -eq 0 ]]; then
    JURIS_LIST=("tas" "act" "vic" "sa" "nt" "wa" "nsw" "qld")
fi

# Commentary 
printf "\nWill run script to organise NGCA files after the AUSPOS processing.\n"
printf "\nWill be working on these jurisdictions: ${JURIS_LIST[*]^^}, for the $ARCHIVE NGCA archive.\n"

# Pause and wait for user reaction
#read -p "Hit [ENTER] key if correct, OR [Ctrl C] if you want to abort..."

for JURIS in ${JURIS_LIST[*]}; do

    cd $SCRIPT_DIR/$JURIS/$ARCHIVE

    # Make backup of files in case there are any issues
    printf "\nCopying files to backup directory in case there are issues later...\n"
    cp -r ../$ARCHIVE/ ../sinexFiles/ ../rinexantls/ ../backup/

    # Coordinate transformation
    # - ITRF2020 --> ITRF2014
    printf "\nTransforming from ITRF2020 to ITRF2014...\n"
    ngca_transformSinex.py -j $JURIS

    # Remove stations outside of GDA2020 from NGCA 
    printf "\nRemoving stations outside GDA2020...\n"
    ngca_exciseStations.py -j $JURIS

    # Select stations from APREF to use as constraints
    printf "\nSelecting stations to use as constraints...\n"
    ngca_selectRef.py -j $JURIS

    # Form baselines from the stations in a SINEX file and create DynaML formatted files
    printf "\nCreating baselines...\n"
    ngca_createBLs.py -j $JURIS

    # Update translation tables in transTables/ with those on the GA ftp server
    printf "\nUpdating transTables...\n"
    ngca_updateTTables.py

    # Check the jurisdicitonal translation tables for duplicates
    printf "\nChecking transTables...\n"
    ngca_checkTTables.py 

    # Translate stations that use the conventional 4-character ID to the name commonly used by the state or jurisdiction
    printf "\nUsing transtables to update station names...\n"
    ngca_stationTrans.py -j $JURIS

    # Run minimally-constrained adjustment for each baseline cluster
    printf "\nAdjusting each cluster\n"
    ngca_minAdjust.py -j $JURIS -c $CPUS

    # Update the Vscale in the NGCA baseline files based on the values in sigma0.dat
    printf "\nUpdating vscale values...\n"
    ngca_addSigma0.py -j $JURIS
    
    # Combine into a single file
    printf "\nCombining clusters...\n"
    ngca_combine.py -j $JURIS -a $ARCHIVE

    # Run fully constrained adjustment performed after NGCA scaling
    printf "\nAdjusting combined NGCA...\n"
    ngca_fullAdjust.py -j $JURIS -a $ARCHIVE

    # Organise results
    printf "\nOrganising NGCA results...\n"
    ngca_results.py -j $JURIS -a $ARCHIVE

done

printf "\n------ DONE ------\n"
printf "\n"

