#!/usr/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$PATH:$SCRIPT_DIR/scripts"

# Archive date for NGCA processing, e.g. the date on which get_ngca.py has downloaded the NGCA archive contents for the jurisdictions being processed 
ARCHIVE=20240715

# List of jurisdictions to run consecutively - increasing number of RINEX
JURIS_LIST=("tas" "act" "vic" "sa" "nt" "wa" "nsw" "qld")

# Alternatively run a single jurisdiction
#JURIS_LIST=("tas")

# Commentary 
printf "\nWill run script to organise NGCA files after the AUSPOS processing.\n"
printf "\nWill be working on these jurisdictions: ${JURIS_LIST[*]^^}, for the $ARCHIVE NGCA archive.\n"

# Pause and wait for user reaction
read -p "Hit [ENTER] key if correct, OR [Ctrl C] if you want to abort..."

for JURIS in ${JURIS_LIST[*]}; do

    # Navigate to jurisdiction directory
    printf "\nNavigating to ${JURIS^^} archive directory...\n"
    cd ~/ngca/$JURIS/$ARCHIVE

    # Make backup of files in case there are any issues
    printf "\nCopying files to backup directory in case there are issues later...\n"
    cp -r ../$ARCHIVE/ ../sinexFiles/ ../rinexantls/ ../backup/

    # renameless 20251111 - Code review has identified script is not required, supporting scripts have been modified to remove namechanges.dat functionality.
    # Revert station names modified by prepfiles
    # printf "\nExecuting nameChanges.pl...\n"
    # nameChanges.pl

    # Coordinate transformation
    # - ITRF2020 --> ITRF2014
    printf "\nExecuting transformSINEX_NGCA.py...\n"
    ngca_transformSinex.py

    # Remove stations outside of GDA2020 from NGCA 
    printf "\nExecuting exciseStationsNGCA.py...\n"
    ngca_exciseStations.py

    # Select stations from APREF to use as constraints
    printf "\nExecuting selectRef.py...\n"
    ngca_selectRef.py

    # Form baselines from the stations in a SINEX file and create DynaML formatted files
    printf "\nExecuting createBLs.py...\n"
    ngca_createBLs.py

    # Update translation tables in transTables/ with those on the GA ftp server
    printf "\nExecuting updateTTables.py...\n"
    ngca_updateTTables.py

    # Check the jurisdicitonal translation tables for duplicates
    printf "\nExecuting checkTTables.py...\n"
    ngca_checkTTables.py 

    # Translate stations that use the conventional 4-character ID to the name commonly used by the state or jurisdiction
    printf "\nExecuting stationTrans.pl...\n"
    ngca_stationTrans.pl

    # Navigate to baseline directory
    cd ../baselines

    # Run minimally-constrained adjustment for each baseline cluster
    printf "\nExecuting minAdjust.pl...\n"
    ngca_minAdjust.pl

    # Report, organise, and generate Sigma0.dat
    printf "\nExecuting howgoesit.pl...\n"
    ngca_howGoesIt.pl

    # Print out the last Sigma0 from the Sigma0 file
    printf "\nHighest Sigma0s for $JURIS:  \n"
    tail -n5 sigma0.dat
    printf "\n"

    # Update the Vscale in the NGCA baseline files based on the values in sigma0.dat
    printf "\nExecuting addSigma0.py...\n"
    ngca_addSigma0.py
    
    # Combine into a single file
    printf "\nExecuting combineNGCA.py...\n"
    ngca_combine.py

    # Run fully constrained adjustment performed after NGCA scaling
    printf "\nExecuting fullAdjust.py...\n"
    ngca_fullAdjust.py

    # Organise results
    printf "\nExecuting ngcaResults.py...\n"
    ngca_results.py

done

printf "\n------ DONE ------\n"
printf "\n"

