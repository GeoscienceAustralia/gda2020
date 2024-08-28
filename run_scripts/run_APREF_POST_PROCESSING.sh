#!/usr/bin/bash

# APREF solution
YYYYMMDD=20240323

# Comments
printf "\nAPREF solution selected: $YYYYMMDD...\n"

# Navigate to APREF working directory
printf "\nNavigating to working directory...\n"
cd /home/fedora/apref/workDir/

# Create symbolic link
ln -s AUS0OPSSNX_*_00U_DSC.SNX soln.snx

# Create file for disconts
cp AUS0OPSSNX_*_00U_DSC.SNX disconts$YYYYMMDD.snx

# Remove stations outside of GDA2020
printf "\nRunning exciseStationsAPREF.py, can take ~2 hour...\n"
exciseStationsAPREF.py

# Transform SINEX coordinates
# - ITRF2020@2015.0 --> ITRF2014@2015.0
# - Comment this out if APREF is ITRF2014
printf "\nRunning transformSINEX_APREF.py, can take some time...\n"
transformSINEX_APREF.py

# Create Type-B uncertanty file
printf "\nRunning createTypeB.py...\n"
createTypeB.py  AUS0OPSSNX_*_00U_SOL.SNX.AUS

# Transform all coordiantes to the GDA2020 epoch
printf "\nRunning sinex2epoch, can take ~6 hours...\n"
sinex2epoch -e20:001 -E1.50379:1.18346:1.20716  AUS0OPSSNX_*_00U_SOL.SNX.AUS &> sinex2epoch.out

# Split APREF into constraint and non-constraint sinex files
printf "\nRunning splitApref.py, can take ~30 minutes...\n"
splitApref.py

# Copy into new filename
cp SNXEPO.SNX.CON.AUS apref$YYYYMMDD.snx

# Create point clusters
printf "\nRunning createPC.pl...\n"
createPC.pl

# Rename disconts file
mv apref*.disconts apref$YYYYMMDD.disconts

# Clean
rm *.xml

# Create baselines
printf "\nRunning aprefCreateBLs.py...\n"
aprefCreateBLs.py $YYYYMMDD

# Reduce size of sinex file
printf "\nRunning reduce_sinex_size.py, can take ~5 minutes...\n"
reduce_sinex_size.py

# Make new APREF solution directory, and move relevant files there
mkdir ../$YYYYMMDD
mv apref* disconts* soln*.snx AUS0OPSSNX_* SNXEPO* typeb.dat sinex2epoch.out ../$YYYYMMDD

# Change up a directory
cd ..

# Run QA script
printf "\nRunning xvsolfin_qa.py...\n"
xvsolfin_qa.py

# Run RVS checking script
printf "\nRunning checkRVS.py...\n"
checkRVS.py

# Move files to APREF solution directory
mv xvsolfin$YYYYMMDD.txt apref$YYYYMMDD.rvs $YYYYMMDD


printf "\n------ DONE ------\n"
printf "\n"
