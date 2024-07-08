#!/usr/bin/bash

# Archive
ARCHIVE=20240515

# Jurisdiction
JURIS="act"

# IP Adress of Geodesy machine (sometimes this changes)
IP_ADRESS=[IP ADRESS]

# Commentary
printf "\nWill run script to copy over files after AUSPOS processing.\n"
printf "\nYou have chosen to work on ${JURIS^^} for the $ARCHIVE NGCA archive.\n"
printf "\n"

# Pause and wait for user reaction
read -p "Hit [ENTER] key if correct, OR [Ctrl C] if you want to abort..."


# Navigate to jurisdiction directory
printf "\nNavigating to ${JURIS^^} NGCA directory...\n"
cd /home/fedora/ngca/$JURIS

# Copy over SINEX files
printf "\nCopying over SINEX files...\n"
scp -r geodesy@${IP_ADRESS}:/data/craig/$ARCHIVE/solutions/*.SNX sinexFiles/

# Copy over RINEX files
printf "\nCopying over RINEX files...\n"
scp -r geodesy@${IP_ADRESS}:/data/craig/$ARCHIVE/solutions/*_ls rinexantls/

# Copy over job submition log file
printf "\nCopying over log_4_submit.txt...\n"
scp -r geodesy@${IP_ADRESS}:/data/craig/$ARCHIVE/log_4_submit.txt $ARCHIVE

# Copy over job dowload file
printf "\nCopying over down_4_SNX.sh...\n"
scp -r geodesy@${IP_ADRESS}:/data/craig/$ARCHIVE/down_4_SNX.sh $ARCHIVE

# Copy over nohup file
printf "\nCopying over nohup.out...\n"
scp -r geodesy@${IP_ADRESS}:/data/craig/$ARCHIVE/nohup.out $ARCHIVE


printf "\n------ DONE ------\n"
printf "\n"

