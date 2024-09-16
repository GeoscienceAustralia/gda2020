#!/usr/bin/env python3

## ===
# DESCRIPTION
# This script is designed to read in an APREF SINEX file 
# that has been cut to the GDA2020 extent.It performs a 
# coordinate transformation to the solution estimates and 
# then will write out a new SINEX file. The input file is 
# copied with .ITRF2020 suffix.
#
# If a TransformationSD object becomes available in GeodePy,
# there is some grayed out code that will handle the transformation
# of uncertainties.
#
# Transformation: ITRF2020@RefEpoch --> ITRF2014@RefEpoch
#
# Input:          AUS0OPSSNX_YYYYDDD_YYYYDDD_00U_SOL.SNX.AUS
#
# Output:         AUS0OPSSNX_YYYYDDD_YYYYDDD_00U_SOL.SNX.AUS 

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
os.chdir(home + '/apref/workDir/')

# Input file
for f in glob.glob('AUS0OPSSNX_*_00U_SOL.SNX.AUS'):
    ifile = f
try:
    ifile
except NameError:
    sys.exit('File not found.')

# Copy original file to new name
os.system('cp ./' + ifile + ' ./' + ifile + '.ITRF2020')

## ===
# READ SINEX BLOCKS
# - Header
# - Comments
# - SITE/ID block
# - SOLUTION/ESTIMATE block

# Header
snx_header = gnss.read_sinex_header_line(ifile)

# File comments
#snx_comments = gnss.read_sinex_custom(ifile, 3, 3)
snx_comments = gnss.read_sinex_comments(ifile)

# SITE/ID Block
snx_block_SiteID = gnss.read_sinex_site_id_block(ifile)

# SOLUTION/EPOCHS Block
snx_block_SolnEpochs = gnss.read_sinex_solution_epochs_block(ifile)

# SOLUTION/ESTIMATE Block
df_SolnEstimate = gnss.sinex2dataframe_solution_estimate(ifile)

# SOLUTION/MATRIX_ESTIMATE Block
# - no TransformationSD object in GeodePy for ITRF2020-->ITRF2014
#matrixEstimate = gnss.read_sinex_matrix(ifile)
#df_vcv = pd.DataFrame(
#                                    matrixEstimate,
#                                    columns=['code', 
#                                            'soln', 
#                                            'xx',
#                                            'xy',
#                                            'yy',
#                                            'xz',
#                                            'yz',
#                                            'zz'],
#)

## ===
# TRANSFORM COORDINATES AND VELOCITIES
# - ITRF2020@2015 --> ITR2014@2015
# - VCV TransformationSD object not available in GeodePy

# Isolate coordinates
Xi = df_SolnEstimate[df_SolnEstimate['par'] == 'STAX']['est'].values
Yi = df_SolnEstimate[df_SolnEstimate['par'] == 'STAY']['est'].values
Zi = df_SolnEstimate[df_SolnEstimate['par'] == 'STAZ']['est'].values

# Isolate velocities
Xvi = df_SolnEstimate[df_SolnEstimate['par'] == 'VELX']['est'].values
Yvi = df_SolnEstimate[df_SolnEstimate['par'] == 'VELY']['est'].values
Zvi = df_SolnEstimate[df_SolnEstimate['par'] == 'VELZ']['est'].values

# Isolate VCVs
# - no TranfomrationSD for list transformation in GeodePy
#VCVi = []
#for i in range(len(df_vcv.code)):
#    
#    # Site VCV
#    Q = np.array([[df_vcv.xx[i], df_vcv.xy[i], df_vcv.xz[i]], 
#                  [df_vcv.xy[i], df_vcv.yy[i], df_vcv.yz[i]],
#                  [df_vcv.xz[i], df_vcv.yz[i], df_vcv.zz[i]]]
#    ) 
#
#    # Append to list of VCVs
#    VCVi.append(Q)

# Transform coordinates and velocities
X = []
Y = []
Z = []
Xv = []
Yv = []
Zv = []
for i in range(len(Xi)):

    # Coordinate transformation
    x, y, z, vcv = transform.conform7(Xi[i], Yi[i], Zi[i], constants.itrf2020_to_itrf2014)

    # Append to list
    X.append(x)
    Y.append(y)
    Z.append(z)

    # Velocity transformation
    xv, yv, zv, vcv = transform.conform7(Xvi[i], Yvi[i], Zvi[i], constants.itrf2020_to_itrf2014_vel)

    # Append to list
    Xv.append(xv)
    Yv.append(yv)
    Zv.append(zv)

## ===
# REPLACE DATAFRAME ESTIMATES
# - Coordinates (X,Y,Z)
# - Velocities (Xv, Yv, Zv)
    
# Coordinates
i = 0
k = 0
j = 0
for r in range(len(df_SolnEstimate.code)):

    # Replace X estimate
    if df_SolnEstimate.loc[r, "par"] == "STAX":
        df_SolnEstimate.loc[r, "est"] = X[i]
        i += 1

    # Replace Y estimate
    if df_SolnEstimate.loc[r, "par"] == "STAY":
        df_SolnEstimate.loc[r, "est"] = Y[k]
        k += 1

    # Replace Z estimate
    if df_SolnEstimate.loc[r, "par"] == "STAZ":
        df_SolnEstimate.loc[r, "est"] = Z[j]
        j += 1

# Velocities
i = 0
k = 0
j = 0
for r in range(len(df_SolnEstimate.code)):

    # Replace X estimate
    if df_SolnEstimate.loc[r, "par"] == "VELX":
        df_SolnEstimate.loc[r, "est"] = Xv[i]
        i += 1

    # Replace Y estimate
    if df_SolnEstimate.loc[r, "par"] == "VELY":
        df_SolnEstimate.loc[r, "est"] = Yv[k]
        k += 1

    # Replace Z estimate
    if df_SolnEstimate.loc[r, "par"] == "VELZ":
        df_SolnEstimate.loc[r, "est"] = Zv[j]
        j += 1

# Write to sinex format
snx_block_SolnEstimate = gnss.dataframe2sinex_solution_estimate(df_SolnEstimate)

## ===
# REPLACE DIAGONAL VCV ELEMENTS OF ORIGINAL VCV with new VCV elements
# 1. Form matrix consisting of individual site VCVs only (from transformed VCVs)
# 2. Create full VCV matrix from original VCV only
# 3. Combine using numpy.where()

# 1) Form matrix of diagonal VCVs only
#Q = gnss.dataframe2matrix_snx_vcv(df_vcv)

# 2)  Create full VCV matrix from original VCV only
# - Form dataframe
# - Form matrix from dataframe
#df_snx_solutionMatrixEstimate = gnss.sinex2dataframe_solution_matrix_estimate(ifile)
#Q0 = gnss.dataframe2matrix_solution_matrix_estimate(df_snx_solutionMatrixEstimate)

# 3) Combine matrices
#    - Replace zero elements with original inter-site-covariances
#    - Convert to lower triangle matrix
#vcv = np.where(Q==0, Q0, Q)
#vcv = np.tril(vcv)

## ===
# WRITE TO SINEX FILE

# Read solution matrix estimate block in SINEX format
snx_block_SolnMatrixEstimate = gnss.read_sinex_solution_matrix_estimate_block(ifile)

# Write SINEX 
gnss.writeSINEX(
                ifile,
                header=snx_header,
                comment=snx_comments, 
                siteID=snx_block_SiteID, 
                solutionEpochs=snx_block_SolnEpochs, 
                solutionEstimate=snx_block_SolnEstimate, 
                solutionMatrixEstimate=snx_block_SolnMatrixEstimate
)
