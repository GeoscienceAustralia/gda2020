## === 
# SETUP

# Packages
using DataFrames
using CSV
using Printf
using Tables
using PrettyTables
using Dates
using EzXML

# Functions
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/dynadjust_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/qaqc_operations.jl")

# MKL paths
ENV["LD_LIBRARY_PATH"] = "/opt/intel/oneapi/mkl/latest/lib:/opt/intel/oneapi/compiler/latest/lib"

# User input
JUR = "ACT"
aprefDate= "20220716"
main_dir = "/home/ec2-user/PROJECTS/GDA2020-QAQC/JADJ"

# Derived variables
jur = lowercase(JUR)
input_dir = "$main_dir/input"
output_dir = "$main_dir/output"

## ===
# ORGANISE JADJ FILES
# - move current files to staging
# - download new files

# Move files to staging
if !isempty(filter(startswith(JUR), readdir(input_dir)))
    @info("$(now()) - MOVE JADJ FILES TO STAGING:                 $JUR")
    files = filter(startswith(JUR), readdir(input_dir))
    run(`mv $input_dir/$files $input_dir/staging/`)
end
if !isempty(filter(startswith(jur), readdir(input_dir)))
    @info("$(now()) - MOVE AUX FILES TO STAGING:                  $JUR")
    files = filter(startswith(jur), readdir(input_dir))
    run(`mv $input_dir/$files $input_dir/staging/`)
end

# Download files
@info("$(now()) - DOWNLOAD FILES:                             $JUR")
run(`aws s3 cp s3://gda2020-ngca/adjustmentData/ $input_dir --recursive --exclude "*" --include "$JUR*"`)
run(`aws s3 cp s3://gda2020-ngca/files/ $input_dir --recursive --exclude "*" --include "$jur*"`)

# Setup jadj station files
@info("$(now()) - SETUP STATION FILE NAMES.")
stn_files = filter(contains(".adj.xml"), readdir(input_dir))
stn_files = filter(contains("GDA2020"), stn_files)
act_stn = filter(startswith("ACT"), stn_files)
nsw_stn = filter(startswith("NSW"), stn_files)
vic_stn = filter(startswith("VIC"), stn_files)
tas_stn = filter(startswith("TAS"), stn_files)
qld_stn = filter(startswith("QLD"), stn_files)
sa_stn = filter(startswith("SA"), stn_files)
nt_stn = filter(startswith("NT"), stn_files)
wa_stn = filter(startswith("WA"), stn_files)
camden_stn = filter(startswith("CAMDEN"), stn_files)
surat_stn = filter(startswith("SURAT"), stn_files)

# Setup jadj measurement files
@info("$(now()) - SETUP MEASUREMENT FILE NAMES.")
msr_files = filter(contains("_msr.xml"), readdir(input_dir))
msr_files = filter(contains("GDA2020"), msr_files)
act_msr = filter(startswith("ACT"), msr_files)
nsw_msr = filter(startswith("NSW"), msr_files)
vic_msr = filter(startswith("VIC"), msr_files)
tas_msr = filter(startswith("TAS"), msr_files)
qld_msr = filter(startswith("QLD"), msr_files)
sa_msr = filter(startswith("SA"), msr_files)
nt_msr = filter(startswith("NT"), msr_files)
wa_msr = filter(startswith("WA"), msr_files)
camden_msr = filter(startswith("CAMDEN"), msr_files)
surat_msr = filter(startswith("SURAT"), msr_files)

## ===
# FORM COMPARISON FILE DATAFRAMES
# i.e. files that contain sites we should exclude from duplicate and near checks. 
# - APREF (constraint)
# - APREF (non-constraint)
# - Ignore
# - Rename
# - Near
 
# APREF (constraint) 
# - From disconts file
df_disconts = disconts2dataframe("$(input_dir)/apref$(aprefDate).disconts")

# APREF non-constraint files
#df_non_constraint = dynstnxml2dataframe("$(input_dir)/apref$(aprefDate)_stn.xml")
#apref_non_constraint_list = df_non_constraint.site
#apref_non_constraint_list = SubString.(apref_non_constraint_list, 1, 4)

# Ignore files
ignore_files = filter(contains(".ignore"), readdir(input_dir))
ignore_list = []
for i in eachindex(ignore_files)
    temp_file = ignore_files[i]
    temp_fp = "$(input_dir)/$(temp_file)"
    jurID = uppercase(split(temp_file, '_')[1])
    if i == 1
        global df_ignore = CSV.read(temp_fp, DataFrame, header=["site"])
        df_ignore = insertcols(df_ignore, 2, :jur=>fill(jurID, length(df_ignore.site)))
    else
        df_temp = CSV.read(temp_fp, DataFrame, header=["site"])
        df_temp = insertcols(df_temp, 2, :jur=>fill(jurID, length(df_temp.site)))
        df_ignore = vcat(df_ignore, df_temp)
    end
end
df_ignore = rstrip.(df_ignore)

# Renaming files
rename_files = filter(contains(".renaming"), readdir(input_dir))
for i in eachindex(rename_files)
    temp_file = rename_files[i]
    temp_fp = "$(input_dir)/$(temp_file)"
    jurID = uppercase(split(temp_file, '_')[1])
    if i == 1
        global df_rename = dynrename2dataframe(temp_fp)
        df_rename = insertcols(df_rename, 3, :jur=>fill(jurID, length(df_rename.rename)))
    else
        df_temp = dynrename2dataframe(temp_fp)
        df_temp = insertcols(df_temp, 3, :jur=>fill(jurID, length(df_temp.rename)))
        df_rename = vcat(df_rename, df_temp)
    end
end
df_rename = rstrip.(df_rename)

# Near files
near_files = filter(endswith(".near"), readdir(input_dir))
for i in eachindex(near_files)
    temp_file = near_files[i]
    temp_fp = "$(input_dir)/$(temp_file)"
    jurID = uppercase(split(temp_file, '_')[1])
    if i == 1
        global df_near = dynnear2dataframe(temp_fp)
        df_near = insertcols(df_near, 3, :jur=>fill(jurID, length(df_near.site)))
    else
        df_temp = dynnear2dataframe(temp_fp)
        df_temp = insertcols(df_temp, 3, :jur=>fill(jurID, length(df_temp.site)))
        df_near = vcat(df_near, df_temp)
    end
end
df_near = rstrip.(df_near)

## ===
# FORM COMBINED JURISDICTION STATION DATAFRAME
# - Used to help identify which jurisdiction to contact
stn_files = filter(contains(".adj.xml"), readdir(input_dir))
for i in eachindex(stn_files)
    temp_file = stn_files[i]
    temp_fp = "$(input_dir)/$(temp_file)"
    jurID = uppercase(split(temp_file, '_')[1])
    if i == 1
        global df_jurStn = dynstnxml2dataframe(temp_fp)
        df_jurStn = insertcols(df_jurStn, 2, :jur=>fill(jurID, length(df_jurStn.site)))
    else
        df_temp = dynstnxml2dataframe(temp_fp)
        df_temp = insertcols(df_temp, 2, :jur=>fill(jurID, length(df_temp.site)))
        df_jurStn = vcat(df_jurStn, df_temp)
    end
end

## ===
# QA/QC STEP 1: DUPLICATE STATION SEARCH
@info("$(now()) - RUN DUPLICATE STATION SEARCH:               $JUR")

# Run DynaAdjust
run(`dnaimport \
    -n gda2020.dup \
    "./input/aprefRename_stn.xml" \
    "./input/aprefRename_msr.xml" \
    "./input/apref$(aprefDate)_stn.xml" \
    "./input/apref$(aprefDate)_msr.xml" \
    "./input/$act_stn" \
    "./input/$nsw_stn" \
    "./input/$vic_stn" \
    "./input/$tas_stn" \
    "./input/$qld_stn" \
    "./input/$sa_stn" \
    "./input/$nt_stn" \
    "./input/$wa_stn" \
    "./input/$act_msr" \
    "./input/$nsw_msr" \
    "./input/$vic_msr" \
    "./input/$tas_msr" \
    "./input/$qld_msr" \
    "./input/$sa_msr" \
    "./input/$nt_msr" \
    "./input/$wa_msr" \
    --prefer-single-x-as-g \
    --flag-unused-stations \
    --ignore-similar-msr \
    --remove-ignored-msr \
    -r GDA2020` 
)

# Read duplicate station file (.dst)
df_dstDup = dyndst2dataframe("$(main_dir)/gda2020.dup.dst")

# Identify genuine duplicates
# - Sites in duplicate list, 
# - but not in discont, ignore, rename lists, or discont[1:4].
# - ToDO: consider the non-constraint file too...
# - ToDo: consider the 4-character code of APREF sites from older APREF.

# Form one list of sites to disregard in duplicate search
non_duplicates = vcat(df_disconts.site, df_ignore.site, df_rename.rename, df_rename.original)
non_duplicates = unique(non_duplicates)

# Genuine duplicates
duplicates = df_dstDup.site[df_dstDup.site .∉ Ref(non_duplicates)]
duplicates = unique(duplicates)

# Form 2nd dataframe with duplicate sites removed due to being an APREF discontinuity
# - i.e. create new dataframe without site names matching CODE_YYYYDOY pattern.
# - There is this check for sites in current APREF, but this will capture 
#   most of the discontinuity sites from jurisdictions using an older APREF.
pattern = r"^[A-Z]{1}[A-Z0-9]{3}_[0-9]{7}\z"
index = []
for i in eachindex(df_dstDup.site)
    name_site = df_dstDup.site[i]
    if occursin(pattern, name_site)
        push!(index, i)
    end
end
df_dstDup1 = df_dstDup[Not(index), :]

# Form one list of sites to disregard in duplicate search
non_duplicates1 = vcat(df_disconts.site, df_ignore.site, df_rename.rename, df_rename.original, SubString.(df_disconts.site, 1, 4))
non_duplicates1 = unique(non_duplicates1)

# Genuine duplicates
duplicates1 = df_dstDup1.site[df_dstDup1.site .∉ Ref(non_duplicates1)]
duplicates1 = unique!(duplicates1)

# Remove files from duplicate check (gda2020.dup*)
dup_files = filter(contains("gda2020.dup"), readdir("./"))
rm.(dup_files)

## ===
# QA/QC STEP 2: NEAR STATION SEARCH
@info("$(now()) - RUN NEAR STATION SEARCH:                    $JUR")

# Run DynaAdjust
run(ignorestatus(`dnaimport \
                -n gda2020.near \
                "./input/aprefRename_stn.xml" \
                "./input/aprefRename_msr.xml" \
                "./input/apref$(aprefDate)_stn.xml" \
                "./input/apref$(aprefDate)_msr.xml" \
                "./input/$act_stn" \
                "./input/$nsw_stn" \
                "./input/$vic_stn" \
                "./input/$tas_stn" \
                "./input/$qld_stn" \
                "./input/$sa_stn" \
                "./input/$nt_stn" \
                "./input/$wa_stn" \
                "./input/$act_msr" \
                "./input/$nsw_msr" \
                "./input/$vic_msr" \
                "./input/$tas_msr" \
                "./input/$qld_msr" \
                "./input/$sa_msr" \
                "./input/$nt_msr" \
                "./input/$wa_msr" \
                --search-nearby-stn \
                --prefer-single-x-as-g \
                --flag-unused-stations \
                --ignore-similar-msr \
                --remove-ignored-msr \
                -r GDA2020`)
)

# Read files
# - Duplicate station file (.dst)
# - Near files (.near)

# DST file 
# - ToDo: extend the dyndst2dataframe to account for near stations too.
df_dst_dup, df_dst_near = dyndst2dataframe("$(main_dir)/gda2020.near.dst")

# Remove files from near check (gda2020.near*)
remove_files = filter(contains("gda2020.near"), readdir("./"))
rm.(remove_files)

# Remove near files due to discontinuities
# - i.e. create new dataframe without site names matching CODE_YYYYDOY pattern
pattern = r"^[A-Z]{1}[A-Z0-9]{3}_[0-9]{7}\z"
index = []
for i in eachindex(df_dst_near.site)

    # Only remove APREF sites
    name_site = df_dst_near.site[i]
    name_nearSite = df_dst_near.nearSite[i]
    if occursin(pattern, name_site) || occursin(pattern, name_nearSite)
        push!(index, i)
    end
end
df_dst_near1 = df_dst_near[Not(index), :]

# Form one list of sites to disregard in near station search
non_near = vcat(df_disconts.site, df_near.site, df_near.nearSite, df_rename.rename, df_rename.original, SubString.(df_disconts.site, 1, 4))
non_near = unique(non_near)

# Genuine near stations
df_nearStns = filter(r -> r.site .∉ Ref(non_near), df_dst_near1)

# Now add jurisdiction columns
# - site jurisdiction. 
# - near site jurisdiction. 
siteJur = []
for site in df_nearStns.site
    j = filter(r->r.site==site, df_jurStn).jur[1]
    push!(siteJur, j)
end
df_nearStns = insertcols(df_nearStns, 2, :siteJur => siteJur)
nearSiteJur = []
for site in df_nearStns.nearSite
    j = filter(r->r.site==site, df_jurStn).jur[1]
    push!(nearSiteJur, j)
end
df_nearStns = insertcols(df_nearStns, 4, :nearSiteJur => nearSiteJur)


## ===
# QA/QC STEP 3: IDENTIFY MARKS WITH BOTH ELLIPSOID AND ORTHOMETRIC HEIGHT
@info("$(now()) - RUN ORTHOMETRIC AND GNSS HEIGHT SEARCH:     $JUR")

# Run DynaAdjust
run(ignorestatus(`dnaimport \
                -n gda2020.m2s \
                "./input/aprefRename_stn.xml" \
                "./input/aprefRename_msr.xml" \
                "./input/apref$(aprefDate)_stn.xml" \
                "./input/apref$(aprefDate)_msr.xml" \
                "./input/$act_stn" \
                "./input/$nsw_stn" \
                "./input/$vic_stn" \
                "./input/$tas_stn" \
                "./input/$qld_stn" \
                "./input/$sa_stn" \
                "./input/$nt_stn" \
                "./input/$wa_stn" \
                "./input/$act_msr" \
                "./input/$nsw_msr" \
                "./input/$vic_msr" \
                "./input/$tas_msr" \
                "./input/$qld_msr" \
                "./input/$sa_msr" \
                "./input/$nt_msr" \
                "./input/$wa_msr" \
                "--output-msr-to-stn" \
                -r GDA2020`)
)

# Read in measurement-to-station file (m2s)
df_m2s, df_m2s_warnings = dynm2s2dataframe("$(main_dir)/gda2020.m2s.m2s")

# Isolate sites into a vector for reference
identified_sites = df_m2s_warnings.site

# Write out new .msr files and edit the ignore flag at identified measurements
for file in msr_files
    @info("$(now()) - WORKING ON M2S WARNINGS:                     $file")
    lines = readlines("$(input_dir)/$(file)")
    f = open("$(output_dir)/$(file)", "w")
    global ignore_flag = false
    for i in eachindex(lines)
        if contains(lines[i], "<Type>H</Type>") && (split(lines[i+4], ['<','>'])[3] in identified_sites || split(lines[i+3], ['<','>'])[3] in identified_sites)
            global ignore_flag = true
            println(f, lines[i])
        elseif ignore_flag == true && contains(lines[i], "<Ignore/>")
            println(f, "    <Ignore>*</Ignore>")
            global ignore_flag = false
        else
            println(f, lines[i])
        end
    end
    close(f)
end
lines = nothing

# Remove files from m2s check (gda2020.m2s*)
remove_files = filter(contains("gda2020.m2s"), readdir("./"))
rm.(remove_files)

## ===
# QA/QC STEP 4: JADJ COORDINATE DELTAS


## ===
# SUMMARY

f = open("$(output_dir)/duplicate_summary_JADJ_QAQC.txt", "w");

@printf(f, "==================\nJADJ QA/QC SUMMARY\n==================\n\n\n")

@printf(f, "--- DUPLICATE STATION SEARCH -------------------------------------------------------------------------\n")
@printf(f, "* Identified by DynAdjust (.DST file).\n")
@printf(f, "* Only genuine duplicates remain.\n\n")
if isempty(duplicates1)

    @printf(f, "    NUMBER OF REMAINING DUPLICATES:     0\n")

else

    @printf(f, "    NUMBER OF REMAINING DUPLICATES:     %i\n\n", length(duplicates1))
    for i in eachindex(duplicates1)
        if i == 1
            global df_duplicates1 = filter(r->r.site==duplicates1[i], df_jurStn)
        else
            df_temp = filter(r->r.site==duplicates1[i], df_jurStn)
            df_duplicates1 = vcat(df_duplicates1, df_temp)
        end
    end
    pretty_table(   f,
                    df_duplicates1,
                    header=(["Site", "Jur", "Constraint", "X", "Y", "Z", "Description"]), 
                    tf=tf_borderless,
    )
end
@printf(f, "\n\n\n")
@printf(f, "--- NEAR STATION SEARCH ------------------------------------------------------------------------------\n")
@printf(f, "* Identified by DynAdjust (.DST file).\n")
@printf(f, "* Only genuine near stations remain.\n\n")
if isempty(df_nearStns.site)

    @printf(f, "    NUMBER OF REMAINING NEAR STATIONS:     0\n")

else
    
    @printf(f, "    NUMBER OF REMAINING NEAR STATIONS:     %i\n\n", length(df_nearStns.site))

    pretty_table(   f,
                    df_nearStns,
                    header=(["Site", "Jur (site)", "Near Site", "Jur (near site)", "ΔHrz", "ΔVrt"]), 
                    tf=tf_borderless,
    )

end
close(f)



print("--- DONE ---")
