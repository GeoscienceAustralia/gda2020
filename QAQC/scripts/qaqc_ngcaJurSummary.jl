## ===
# DESCRIPTION
# - This script will collate relevant details to summarise an NGCA archive. It is designed to be run
#   after the NGCA processing is complete. Its purpose is to aid the troubleshooting of NCGA processing.
# - Unless archive already exists in the 'archives' directory, the NGCA solution will be downloaded
#   from s3://gda2020-ngca/ngcaResults/.
# 
# Details summarised:
#                       - Size of archive.
#                       - Number of RINEX.
#                       - Number of Clusters.
#                       - Number of clusters not processed.
#                       - Number of SINEX.
#                       - DynAdjust result summary and ten largest baseline corrections.
#                       - Ten largest coordinate uncertainties from AUSPOS solutions.
#                       - Ten largest coordinate differences between archives.
#                       - Ten largest coordinate uncertainty differences between archives.
#
# Input:
#                       - Two NGCA solutions from the same jurisdiction (.zip).
#
# Output:
#                       - Text file of summary (JUR_NGCA_YYYYMMDD_SUMMARY.txt).
#                       - CSV file of collected SINEX data (JUR_NGCA_YYYYMMDD_SUMMARY.csv).
#                       - DAT file of cluster files (JUR_NGCA_YYYYMMDD_CLUSTERS.dat).
#                       - Optional: Text file of all coodinate differences between archives 
#                         (JUR_NGCA_COORDINATE_DIFFERENCES_YYYYMMDD_YYYYMMDD.txt). 


## === 
# SETUP

# Packages
using DataFrames
using CSV
using StatsBase
using Printf
using Tables
using PrettyTables
using Dates

# Functions
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/sinex_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/dynadjust_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/qaqc_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/geodetic_operations.jl")

# User input
# - Jurisdiction
# - Archive date
# - Main directory path
# - Archive comparison condition
# - Save archive comparison file
# - Dictionary of jurisdiction representatives
JUR = "NSW"
date = "20250115"
date2 = "20241115"
main_dir = "/home/ec2-user/PROJECTS/GDA2020-QAQC/NGCA/PROD"
archive_compare = true
saveCoordinateDeltaFile = true
JUR_REP = Dict(["ACT"=>"Gavin", 
                "TAS"=>"Matthew", 
                "VIC"=>"Josh", 
                "SA"=>"Marc", 
                "NT"=>"Bill", 
                "WA"=>"Khandu", 
                "QLD"=>"Darren",
                "NSW"=>"Joel"]
)

# Derived variables
# - archives directory path
# - archive folder name
# - output path
# - List of jurisdiction representatives
jur = lowercase(JUR)
archive_name = "$(JUR)_NGCA_$(date)"
archive_dir = "$(main_dir)/archives/$(archive_name)"
archive_name2 = "$(JUR)_NGCA_$(date2)"
archive_dir2 = "$(main_dir)/archives/$(archive_name2)"
output_dir = "$(main_dir)/output"
s3_dir = "s3://gda2020-ngca/ngcaResults"

@info("$(now()) - === NGCA SUMMARY BEGIN ===")

# Download archive
if isfile("$(archive_dir).zip")
    @info("$(now()) - NGCA ZIP FOLDER:                            $(archive_name).zip | already exists\n")
else
    @info("$(now()) - NGCA ZIP FOLDER:                            $(archive_name).zip | downloading\n")
    run(`aws s3 cp $(s3_dir)/$(archive_name).zip $(main_dir)/archives/`)
end
if isdir("$(archive_dir)")
    @info("$(now()) - NGCA UNZIPPED FOLDER:                       $(archive_name) | already exists")
else
    @info("$(now()) - NGCA UNZIPPED FOLDER:                       $(archive_name) | extracting\n")
    run(`unzip $(archive_dir).zip -d $(main_dir)/archives/`)
end

@info("$(now()) - INITIAL ARCHIVE QUERIES.")

## ===
# ARCHIVE SIZE
# - MegaBytes
zip_size = stat("$(archive_dir).zip").size/1000000

## ===
# CLUSTERS
# - List of cluster names
# - Number of clusters
cluster_names = filter(endswith("_ls"), readdir("$archive_dir/clusters"))
num_clusters = length(cluster_names)

## ===
# RINEX
# - List of rinex from clusters (rnx_list)
# - Number of rinex from clusters (num_rnx)
# - Count rinex occurences (rnx_count)
# - Dataframe of rinex count (df_rnx_count)
# - Number of rinex names that occur more than once (num_rnx_duplicates)
# - Get mark names for use later
rnx_list = String[]
for f in cluster_names
    rnx = CSV.read("$archive_dir/clusters/$f", DataFrame, delim=' ', header=false, ignorerepeated=false)[:,1]  
    global rnx_list = vcat(rnx_list, rnx)
end
num_rnx = length(rnx_list)
rnx_count = countmap(rnx_list)
k = collect(keys(rnx_count))
v = collect(values(rnx_count))
df_rnx_count = DataFrame(rinex=k, count=v)
df_rnx_duplicates = filter(r -> r.count > 1, df_rnx_count)
num_rnx_duplicates = length(df_rnx_duplicates.rinex)
marker_names = SubString.(rnx_list, 1,4)
marker_names = unique(marker_names)

## ===
# SINEX
# - List of SINEX names
# - Number of SINEX
snx_names = filter(endswith(".SNX.$(JUR).NGCA"), readdir("$archive_dir/snx"))
num_snx = length(snx_names)

## ===
# CLUSTERS NOT PROCESSED
# - Get list of sinex codes
# - Get list of cluster codes
# - Identify which cluster codes are not in sinex code list
# - Identify AUSPOS job numbers associated with those cluster codes

@info("$(now()) - IDENTIFYING CLUSTERS NOT PROCESSED.")

if num_clusters-num_snx != 0
    snx_code = []
    for i in eachindex(snx_names)
        code = split(snx_names[i], ".")[1]
        global snx_code = push!(snx_code, code)
    end
    cluster_code = []
    for i in eachindex(cluster_names)
        # Comment in/out for single mark
        code = split(cluster_names[i], "_")[1]
        #code = first(cluster_names[i], 11) 
        global cluster_code = push!(cluster_code, code)
    end
    clusters_not_processed = cluster_code[cluster_code .∉ Ref(snx_code)]
    job_numbers_not_processed = []
    for i in eachindex(clusters_not_processed)
        job_lines = readlines("$(archive_dir)/feedback/log_4_submit.txt")
        job_line_index = findfirst(contains(" $(clusters_not_processed[i])_ls "), job_lines)
        job_line = job_lines[job_line_index]
        job = job_line[1:22]
        push!(job_numbers_not_processed, job)
    end
end

## ===
# RINEX FLAGGED
# - Flagged by verifySub.py

@info("$(now()) - IDENTIFYING RINEX FLAGGED.")

if isfile("$(archive_dir)/feedback/verifySub.log")
    rnx_flagged = CSV.read("$archive_dir/feedback/verifySub.log", DataFrame, delim=':', header=false, ignorerepeated=true)[:,3]
    num_rnx_flagged = length(rnx_flagged)
    verifySub_strings = CSV.read("$archive_dir/feedback/verifySub.log", DataFrame, delim=':', header=false, ignorerepeated=true)[:,4]
else
    num_rnx_flagged = 0
end

## ===
# RINEX NAME CHANGES
# - create dictionary of name changes
# - only if nameChanges.dat exists
if isfile("$archive_dir/other/nameChanges.dat")
    @info("$(now()) - IDENTIFYING RINEX NAME CHANGES.")
    df_name_changes = CSV.read("$archive_dir/other/nameChanges.dat", DataFrame, delim=' ', header=false, ignorerepeated=false)
    rnx_name_changes = Dict(Pair.(df_name_changes[:,1], df_name_changes[:,2]))
end

## ===
# ADJUSTMENT RESULTS 
# - From JUR_GDA2020.phased.adj
# - Adjustment summary block
# - Adjusted coordinates

@info("$(now()) - READ ADJUSTMENT FILE:                       $(JUR)_GDA2020.phased.adj")

# Read adjustment summary block
adjustment_summary_block = read_DynAdjustAdjustmentSummary_BLOCK("$(archive_dir)/adjustments/$(JUR)_GDA2020.phased.adj")

# Adjusted measurements
# - Filter to baselines only (X measurement type)
# - Form new baseline based measurement correction dataframe
# - Add 3D correction
df_adjusted_msr = dynadj2dataframe("$(archive_dir)/adjustments/$(JUR)_GDA2020.phased.adj")

# Baseline measurements only
df_adjusted_msr = filter(r -> r.type=="X", df_adjusted_msr)

# Isolate X,Y,Z corrections
dfx = filter(r -> r.par == "X", df_adjusted_msr)
dfy = filter(r -> r.par == "Y", df_adjusted_msr)
dfz = filter(r -> r.par == "Z", df_adjusted_msr)

# Form new dataframe
df_adjusted_msr = DataFrame(
                            station1=dfx.station1,
                            station2=dfx.station2,
                            corX=dfx.cor,
                            corY=dfy.cor,
                            corZ=dfz.cor,
                            nstatX= dfx.nstat,
                            nstatY= dfy.nstat,
                            nstatZ= dfz.nstat,
)
df_adjusted_msr = insertcols(df_adjusted_msr, :cor3D => sqrt.(df_adjusted_msr.corX.^2 + df_adjusted_msr.corY.^2 + df_adjusted_msr.corZ.^2))

# Clean
dfx = nothing
dfy = nothing
dfz = nothing

## ===
# AUSPOS COORDINATE SUMMARY
# - Read each SINEX into a dataframe
# - Combine into one dataframe
# - Remove non-jurisdiction sites
# - Add marker number

@info("$(now()) - AUSPOS COORDINATE SUMMARY.")

df = sinex2dataframe_NGCA(archive_dir, JUR, marker_names)

# Save file
output_file = "$(JUR)_NGCA_$(date)_SUMMARY"
CSV.write("$output_dir/$(output_file).csv", df)

# ===
# SUBSET DATAFRAME
# - For display purposes
df_sigmas = select(df, :site, :number, :meanEpoch, :sigmaE, :sigmaN, :sigmaU)
df_sigmas = insertcols(df_sigmas, :sigma3D => sqrt.(df_sigmas.sigmaE.^2 + df_sigmas.sigmaN.^2 + df_sigmas.sigmaU.^2))

## ===
# CLUSTER DATAFRAME
# - For saving to file
df_clusters = DataFrame([[],[],[],[]], ["rinex","antenna","antHgt","cluster"])
for i in eachindex(cluster_names)
    cluster_name = cluster_names[i]
    local df_temp = CSV.read("$archive_dir/clusters/$cluster_name", DataFrame, delim=' ', header=[:rinex, :antenna, :antHgt], ignorerepeated=false)
    cluster_col = fill(cluster_name, length(df_temp[:,1]))
    df_temp = insertcols(df_temp, :cluster => cluster_col)
    global df_clusters = vcat(df_clusters, df_temp)
end
df_clusters = select(df_clusters, :cluster, :rinex, :antenna, :antHgt)
# Add job numbers
job_numbers = []
for i in eachindex(df_clusters.cluster)
    job_lines = readlines("$(archive_dir)/feedback/log_4_submit.txt")
    job_line_index = findfirst(contains(" $(df_clusters.cluster[i]) "), job_lines)
    job_line = job_lines[job_line_index]
    job = job_line[1:22]
    push!(job_numbers, job)
end
df_clusters.job = job_numbers
# Save to file
output_file = "$(JUR)_NGCA_$(date)_CLUSTERS"
fid = open("$output_dir/$output_file.dat","w")
for i in eachindex(df_clusters.cluster)
    write(fid, @sprintf("%9s  %12s  %14s  %.4f  %s\n",  df_clusters.cluster[i], df_clusters.rinex[i], df_clusters.antenna[i], df_clusters.antHgt[i], df_clusters.job[i]))
end
close(fid)

## ===
# ARCHIVE COORDINATE COMPARISON (optional) 
# - Read each SINEX into a dataframe
# - Combine into one dataframe
# - Remove non-jurisdiction sites
# - Add marker number
# - Compute coordinate deltas
# - Compute uncertainty deltas

if archive_compare == true

    @info("$(now()) - ARCHIVE COORDINATE COMPARISON.")

    # Archive size
    zip_size2 = stat("$(archive_dir2).zip").size/1000000

    # Setup
    cluster_names2 = filter(endswith("_ls"), readdir("$archive_dir2/clusters"))
    rnx_list2 = String[]
    for f in cluster_names2
        rnx2 = CSV.read("$archive_dir2/clusters/$f", DataFrame, delim=' ', header=false, ignorerepeated=false)[:,1]  
        global rnx_list2 = vcat(rnx_list2, rnx2)
    end
    marker_names2 = SubString.(rnx_list2, 1,4)
    marker_names2 = unique(marker_names2)

    df2 = sinex2dataframe_NGCA(archive_dir2, JUR, marker_names2)

    # Coordinate deltas
    df_archiveDelta = coordinate_delta_SINEX(df, df2)
    
    # Uncertainty deltas
    df_archiveUncertaintyDelta = uncertainty_delta_SINEX(df, df2)

end

## === 
# DISPLAY ARCHIVE SUMMARY 
# - Email preamble
# - Archive summary
# - Adjustment summary
# - AUSPOS coordinate uncertainty summary

@info("$(now()) - WRITE ARCHIVE SUMMARY TO FILE:              $(JUR)_NGCA_$(date)_SUMMARY.txt")

f = open("$(output_dir)/$(JUR)_NGCA_$(date)_SUMMARY.txt", "w");

@printf(f, "\nHello %s, \
            \n\nThe above NGCA solution is now available via SFTP on the S3 bucket. \
            \n\nLet us know if there are issues. \
            \n\nBelow is summary information. It's only there to aid troubleshooting if required. \
            \n\nThanks, \
            \nGDA2020 Team\n\n", JUR_REP["$(JUR)"])

@printf(f, "\n--- SUMMARY ------------------------------------------------------------------------------------------\n\n")

@printf(f, "Archive name:                                      %s\n", archive_name)
if archive_compare == true
    @printf(f, "Archive size:                                      %.2f MB (compared to %.2f MB of %s)\n", zip_size, zip_size2, archive_name2)
else
    @printf(f, "Archive size:                                      %.2f MB\n", zip_size)
end
@printf(f, "Number of valid RINEX:                             %i\n", num_rnx)
@printf(f, "Number of duplicate RINEX:                         %i\n", num_rnx_duplicates)
@printf(f, "Number of RINEX flagged by verifySub.py:           %i\n", num_rnx_flagged)
@printf(f, "Number of Clusters:                                %i\n", num_clusters)
@printf(f, "Number of SINEX:                                   %i\n", num_snx)
@printf(f, "Number of clusters not processed:                  %i\n", num_clusters - num_snx)

if num_clusters-num_snx != 0
    @printf(f, "\n--- UNPROCESSED CLUSTERS -----------------------------------------------------------------------------\n\n")
    for i in eachindex(clusters_not_processed)
        @printf(f, "   %9s_ls - %s (AUSPOS job number)\n", clusters_not_processed[i], job_numbers_not_processed[i])
        cluster_string = clusters_not_processed[i]
        cluster_path = "$(archive_dir)/clusters/$(cluster_string)_ls"
        cluster_lines = readlines(cluster_path)
        if isfile("$archive_dir/other/nameChanges.dat")
            for j in eachindex(cluster_lines)
                rnx_string = SubString(cluster_lines[j], 1,12)
                if haskey(rnx_name_changes, rnx_string)
                    @printf(f, "                  %s - (%s)\n", cluster_lines[j], rnx_name_changes[rnx_string])
                else
                    @printf(f, "                  %s\n", cluster_lines[j])
                end
            end
        else
            for j in eachindex(cluster_lines)
                @printf(f, "                  %s\n", cluster_lines[j])
            end
        end
    end
end

if num_rnx_duplicates > 0
    @printf(f, "\n--- DUPLICATES ---------------------------------------------------------------------------------------\n\n")
    pretty_table(df_rnx_duplicates)
end
if num_rnx_flagged != 0
    @printf(f, "\n--- RINEX FLAGGED ------------------------------------------------------------------------------------\n")
    @printf(f, "* From verifySub.log\n\n")
    for i in eachindex(rnx_flagged)

        if isfile("$archive_dir/other/nameChanges.dat")
            if haskey(rnx_name_changes, rnx_flagged[i])
                @printf(f, "   %s -%s - (%s)\n", rnx_flagged[i], verifySub_strings[i], rnx_name_changes[rnx_flagged[i]])
            else
                @printf(f, "   %s -%s\n", rnx_flagged[i], verifySub_strings[i])
            end
        else
            @printf(f, "   %s -%s\n", rnx_flagged[i], verifySub_strings[i])
        end
    end
end

@printf(f, "\n\n--- COMBINED ADJUSTMENT RESULTS ----------------------------------------------------------------------\n")
@printf(f, "* From DynAdjust solution\n\n")

for i in eachindex(adjustment_summary_block)
    @printf(f, "%s\n", adjustment_summary_block[i])
end

@printf(f, "\n\nTen largest baseline corrections:\n")
@printf(f, "* Units: metres\n\n")
pretty_table(
                f, 
                sort(df_adjusted_msr, order(:cor3D, rev=true))[1:10,:], 
                show_subheader=false, 
                formatters=ft_printf("%.5f"), 
                tf=tf_borderless,
)

@printf(f, "\n--- TEN LARGEST COORDINATE UNCERTAINTIES -------------------------------------------------------------\n")
@printf(f, "* From AUSPOS solution\n")

pretty_table(
                f, 
                sort(df_sigmas, order(:sigma3D, rev=true))[1:10,:],
                header=(["Site", "Number", "Date", "σE  ", "σN  ", "σU  ", "σ3D "], 
                ["(code)", "(marker)", "(mean epoch)", "(m)  ", "(m)  ", "(m)  ", "(m)  "]),
                formatters=ft_printf("%.5f"), 
                tf=tf_borderless,
                alignment=[:l, :l, :c, :r, :r, :r, :r],
)

if archive_compare == true

    @printf(f, "\n--- TEN LARGEST COORDINATE DIFFERENCES ---------------------------------------------------------------\n")
    @printf(f, "* From AUSPOS solution\n")
    @printf(f, "* In comparison to %s\n\n", "$(JUR)_NGCA_$(date2)")

    pretty_table(
                    f, 
                    sort(select(df_archiveDelta, :site, :soln, :dE, :dN, :dU, :distance3D), 
                    order(:distance3D, rev=true))[1:10,:],
                    header=(["Site", "Soln", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                    [" ", " ", "(mm)", "(mm)", "(mm)", "(mm)"]),
                    formatters=ft_printf("%.2f"), 
                    tf=tf_borderless,
                    alignment=[:l, :c, :r, :r, :r, :r],
    )

    @printf(f, "\n--- TEN LARGEST UNCERTAINTY DIFFERENCES --------------------------------------------------------------\n")
    @printf(f, "* From AUSPOS solution\n")
    @printf(f, "* In comparison to %s\n\n", "$(JUR)_NGCA_$(date2)")

    pretty_table(
                    f, 
                    sort(select(df_archiveUncertaintyDelta, :site, :soln, :dE, :dN, :dU, :distance3D), 
                    order(:distance3D, rev=true))[1:10,:],
                    header=(["Site", "Soln", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                    [" ", " ", "(mm)", "(mm)", "(mm)", "(mm)"]),
                    formatters=ft_printf("%.2f"), 
                    tf=tf_borderless,
                    alignment=[:l, :c, :r, :r, :r, :r],
    )

end

@printf(f, "\n\n------------------------------------------------------------------------------------------------------\n\n")

close(f)

if saveCoordinateDeltaFile == true && archive_compare == true

    # Add in job number column
    job_numbers = []
    for i in eachindex(df_archiveDelta.cluster)
        cluster = df_archiveDelta.cluster[i]
        job = filter(r->r.cluster==cluster, df_clusters)[1, "job"]
        push!(job_numbers, job)
    end
    df_archiveDelta.job = job_numbers

    f = open("$(output_dir)/$(JUR)_NGCA_COORDINATE_DIFFERENCES_$(date)_$(date2).txt", "w");
    pretty_table(
                        f, 
                        sort(select(df_archiveDelta, :site, :soln, :cluster, :job, :dE, :dN, :dU, :distance3D), 
                        order(:distance3D, rev=true)),
                        header=(["site", "soln", "cluster", "job", "ΔE", "ΔN", "ΔU", "Δ3D"], 
                        [" ", " ", " ", " ", "(mm)", "(mm)", "(mm)", "(mm)"]),
                        formatters=ft_printf("%.2f"), 
                        tf=tf_borderless,
                        alignment=[:l, :l, :l, :l, :r, :r, :r, :r],
    )
    close(f)
end

@info("$(now()) - === NGCA SUMMARY COMPLETE ===")