## ===
# DESCRIPTION
# This script will download a single IGS cumulative solution 
# file from CDDIS. The existance of the file in the local directory 
# will be checked first. See CDDIS download information for access
# setup. The SINEX file will also be converted to a Dataframe, and
# then saved to CSV. 
#
# Prerequisite:    .netrc file in ~/ 


## ===
# SETUP

# Packages
using Dates
using DataFrames
using Tables
using CSV
using Printf

# Functions
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/sinex_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/geodetic_operations.jl")

# User input
gps_week = "2309"
igs_file = "IGS0OPSSNX_1994002_2024104_00U_CRD.SNX"
dsc_file = "IGS0OPSSNX_1994002_2024104_00U_DSC.SNX"
main_dir = "/home/ec2-user/PROJECTS/GDA2020-QAQC/APREF/PROD"

# Derived variables
igs_dir = "$(main_dir)/igs"
igs_fp = "$(igs_dir)/$(igs_file)"
dsc_fp = "$(igs_dir)/$(dsc_file)"

## ===
# DOWNLOAD
# - Confirm file does not already exist
# - Then download and uncompress

# SOLUTION
if isfile("$igs_fp")
    @info("$(now()) - IGS SINEX FILE:                           $igs_file | already exists")
else
    @info("$(now()) - IGS SINEX FILE:                           $igs_file | downloading and extracting")
    run(`wget --auth-no-challenge https://cddis.nasa.gov/archive/gnss/products/$gps_week/$igs_file.Z -P $igs_dir/`)
    run(`gunzip $igs_dir/$igs_file.Z`)
end

# DISCONTS
if isfile("$dsc_fp")
    @info("$(now()) - IGS SINEX FILE:                           $dsc_file | already exists")
else
    @info("$(now()) - IGS SINEX FILE:                           $dsc_file | downloading and extracting")
    run(`wget --auth-no-challenge https://cddis.nasa.gov/archive/gnss/products/$gps_week/$disconts_file.Z -P $igs_dir/`)
    run(`gunzip $igs_dir/$dsc_file.Z`)
end

## ===
# SAVE DATAFRAME AS CSV
# - Using sinex2dataframe()

sinex2dataframe(igs_fp, igs_dir, false, dsc_fp)