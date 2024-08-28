#!/usr/bin/env python3

## ===
# DESCRIPTION
# This script is designed to read in all of the SINEX files
# from the latest NGCA archive of a particular jurisdiction
# It performs a coordinate transformation to the solution 
# estimates and then will write out a new SINEX file. The 
# input file is copied with .ITRF2020 suffix.
#
# Transformation: ITRF2020@ObsEpoch --> ITRF2014@ObsEpoch
#
# Input:          [sinex_name].SNX
#
# Output:         [sinex_name].SNX
#
# To Do:
#           - If a TransformationSD object becomes available in GeodePy, there 
#             is some grayed out code in transformSINEX_APREF.py which could 
#             be used.

## ===
# SETUP

# Packages
import glob
import os
import sys
import pandas as pd
import numpy as np
from geodepy import gnss, transform, constants

# Directory
# - Get the path to the home directory and move to the working directory
home = os.path.expanduser("~")
os.chdir(home + '/PROJECTS/SINEX/NGCA_AUSPOS_ITRF2014_TO_ITRF2020/working')

# Input file
for f in glob.glob('*.SNX'):

    ifile = f

    # Copy original file to new name
    new_name = ifile.split('.')[0]
    os.system('cp ./' + ifile + ' ./' + new_name + '.ITRF2020')

    ## ===
    # READ SINEX BLOCKS
    # - Header
    # - FILE/REFERENCE
    # - INPUT/ACKNOWLEDGMENTS
    # - SOLUTION/STATISTICS
    # - SITE/ID
    # - SITE/RECEIVER
    # - SITE/ANTENNA
    # - SITE/GPS_PHASE_CENTER
    # - SITE/ECCENTRICITY
    # - SOLUTION/EPOCHS
    # - SOLUTION/ESTIMATE (dataframe for computation)
    # - SOLUTION/APRIORI (dataframe for computation)
    # - SOLUTION/MATRIX_ESTIMATE
    # - SOLUTION/MATRIX_APRIORI

    # Header
    snx_header = gnss.read_sinex_header_line(ifile)

    # FILE/REFERENCE
    snx_fileReference = gnss.read_sinex_file_reference_block(ifile)

    # INPUT/ACKNOWLEDGMENTS
    snx_inputAcknowledgements = gnss.read_sinex_input_acknowledgments_block(ifile)
    
    # SOLUTION/STATISTICS
    snx_solnStatistics = gnss.read_sinex_solution_statistics_block(ifile)

    # SITE/ID
    snx_siteID = gnss.read_sinex_site_id_block(ifile)

    # SITE/RECEIVER
    snx_siteReceiver = gnss.read_sinex_site_receiver_block(ifile)

    # SITE/ANTENNA
    snx_siteAntenna = gnss.read_sinex_site_antenna_block(ifile)

    # SITE/GPS_PHASE_CENTER
    snx_siteGpsPhaseCenter = gnss.read_sinex_site_gps_phase_center_block(ifile)

    # SITE/ECCENTRICITY
    snx_siteEccentricity = gnss.read_sinex_site_eccentricity_block(ifile)

    # SOLUTION/EPOCHS
    snx_solnEpochs = gnss.read_sinex_solution_epochs_block(ifile)

    # SOLUTION/ESTIMATE 
    # - dataframe
    df_solnEstimate = gnss.sinex2dataframe_solution_estimate(ifile)

    # SOLUTION/APRIORI
    # - dataframe
    df_solnApriori = gnss.sinex2dataframe_solution_apriori(ifile)

    # SOLUTION/MATRIX_ESTIMATE
    snx_block_solnMatrixEstimate = gnss.read_sinex_solution_matrix_estimate_block(ifile)

    # SOLUTION/MATRIX_APRIORI
    snx_block_solnMatrixApriori = gnss.read_sinex_solution_matrix_apriori_block(ifile)

    ## ===
    # TRANSFORM COORDINATES (Solution Estimate)
    # - ITRF2020@ObsEpoch --> ITR2014@ObsEpoch

    # Isolate coordinates
    Xi = df_solnEstimate[df_solnEstimate['par'] == 'STAX']['est'].values
    Yi = df_solnEstimate[df_solnEstimate['par'] == 'STAY']['est'].values
    Zi = df_solnEstimate[df_solnEstimate['par'] == 'STAZ']['est'].values

    # Transform coordinates
    X = []
    Y = []
    Z = []
    for i in range(len(Xi)):
    
        # Coordinate transformation
        x, y, z, vcv = transform.conform7(Xi[i], Yi[i], Zi[i], constants.itrf2020_to_itrf2014)
    
        # Append to list
        X.append(x)
        Y.append(y)
        Z.append(z)

    ## ===
    # REPLACE DATAFRAME ESTIMATES (Solution Estimate)
    # - Coordinates (X,Y,Z)
    
    # Coordinates
    i = 0
    k = 0
    j = 0
    for r in range(len(df_solnEstimate.code)):

        # Replace X estimate
        if df_solnEstimate.loc[r, "par"] == "STAX":
            df_solnEstimate.loc[r, "est"] = X[i]
            i += 1
    
        # Replace Y estimate
        if df_solnEstimate.loc[r, "par"] == "STAY":
            df_solnEstimate.loc[r, "est"] = Y[k]
            k += 1
    
        # Replace Z estimate
        if df_solnEstimate.loc[r, "par"] == "STAZ":
            df_solnEstimate.loc[r, "est"] = Z[j]
            j += 1

    # Write to sinex format
    snx_block_solnEstimate = gnss.dataframe2sinex_solution_estimate(df_solnEstimate)

    ## ===
    # TRANSFORM COORDINATES (Solution Apriori)
    # - ITRF2020@ObsEpoch --> ITR2014@ObsEpoch

    # Isolate coordinates
    Xi = df_solnApriori[df_solnApriori['par'] == 'STAX']['est'].values
    Yi = df_solnApriori[df_solnApriori['par'] == 'STAY']['est'].values
    Zi = df_solnApriori[df_solnApriori['par'] == 'STAZ']['est'].values

    # Transform coordinates
    X = []
    Y = []
    Z = []
    for i in range(len(Xi)):
    
        # Coordinate transformation
        x, y, z, vcv = transform.conform7(Xi[i], Yi[i], Zi[i], constants.itrf2020_to_itrf2014)
    
        # Append to list
        X.append(x)
        Y.append(y)
        Z.append(z)

    ## ===
    # REPLACE DATAFRAME ESTIMATES (Solution Apriori)
    # - Coordinates (X,Y,Z)
    
    # Coordinates
    i = 0
    k = 0
    j = 0
    for r in range(len(df_solnApriori.code)):

        # Replace X estimate
        if df_solnApriori.loc[r, "par"] == "STAX":
            df_solnApriori.loc[r, "est"] = X[i]
            i += 1
    
        # Replace Y estimate
        if df_solnApriori.loc[r, "par"] == "STAY":
            df_solnApriori.loc[r, "est"] = Y[k]
            k += 1
    
        # Replace Z estimate
        if df_solnApriori.loc[r, "par"] == "STAZ":
            df_solnApriori.loc[r, "est"] = Z[j]
            j += 1

    # Write to sinex format
    snx_block_solnApriori = gnss.dataframe2sinex_solution_apriori(df_solnApriori)

    ## ===
    # WRITE TO SINEX FILE

    # Write SINEX 
    gnss.writeSINEX(
                    ifile,
                    header=snx_header,
                    fileReference=snx_fileReference, 
                    inputAcknowledgments=snx_inputAcknowledgements, 
                    solutionStatistics=snx_solnStatistics,
                    siteID=snx_siteID,
                    siteReceiver=snx_siteReceiver,
                    siteAntenna=snx_siteAntenna,
                    siteGpsPhaseCenter=snx_siteGpsPhaseCenter,
                    siteEccentricity=snx_siteEccentricity,
                    solutionEpochs=snx_solnEpochs,
                    solutionEstimate=snx_block_solnEstimate,
                    solutionApriori=snx_block_solnApriori, 
                    solutionMatrixEstimate=snx_block_solnMatrixEstimate,
                    solutionMatrixApriori=snx_block_solnMatrixApriori,
    )
