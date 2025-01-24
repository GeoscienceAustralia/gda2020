## ===
# DESCRIPTION
# - This script will summarise a national adjustment solution (NADJ). 
# - It is designed to be run after the NADJ processing is complete. 
# - Its purpose is to monitor any changes and aid troubleshooting. 
# 
# Details summarised:   - Processing details.
#                       - Input files.
#                       - Adjustment statistics.
#                       - Ten largest coordinate corrections.
#                       - Largest n-statistics.
#                       - Ten largest coordinate changes (NADJ).
#                       - Ten largest uncertainty changes (NADJ).
#                       - Ten largest coordinate changes (RVS).
#                       - Stations added/removed.
#
# Input:                - NADJ solution 1 (gda2020_YYYYDMMDD.zip).
#                       - NADJ solution 2 (gda2020_YYYYDMMDD.zip).
#                       - RVS file (RVS_GDA2020.txt).
#                       
#
# Output:               - Summary text file (GDA2020_NADJ_YYYYMMDD_YYYYMMDD_SUMMARY.txt).
#                       - Probability density function of cordinate differences (.png).
#                       - Probability density function of uncertainty differences (.png).
#                       - Bar graph of coordinate differences (.png). 

## === 
# SETUP

# Packages
using DataFrames
using CSV
using Printf
using Tables
using PrettyTables
using Dates
using StatsPlots
using Statistics

# Functions
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/geodetic_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/qaqc_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/dynadjust_operations.jl")

# User input
# - date of nadj solution 1
# - date of nadj solution 2
# - working directory (main)
date1 = "20241202"
date2 = "20241001"
rvs_file = "RVS_GDA2020.txt"
main_dir = "/home/ec2-user/PROJECTS/GDA2020-QAQC/NADJ/PROD"

# Derived variables
nadj_archive1 = "gda2020_$(date1)"
nadj_archive2 = "gda2020_$(date2)"
input_dir = "$(main_dir)/input"
output_dir = "$(main_dir)/output"
nadj_dir1 = "$(input_dir)/$(nadj_archive1)"
nadj_dir2 = "$(input_dir)/$(nadj_archive2)"
xyz_fp1 = "$(nadj_dir1)/$(nadj_archive1).phased-stage.xyz"
xyz_fp2 = "$(nadj_dir2)/$(nadj_archive2).phased-stage.xyz"
rvs_fp = "$(main_dir)/other/$(rvs_file)"

@info("$(now()) - === NADJ SUMMARY BEGIN ===")
@info("$(now()) - INPUT NADJ FOLDER (1):                      $nadj_archive1")
@info("$(now()) - INPUT NADJ FOLDER (2):                      $nadj_archive2")
@info("$(now()) - RVS FILE:                                   $rvs_file")

## ===
# DOWNLOAD SOLUTIONS

# NADJ solution (1)
if isfile("$(nadj_dir1).zip")
    @info("$(now()) - NADJ ZIP FOLDER (solution 1):               $(nadj_archive1).zip | already exists")
else
    @info("$(now()) - NADJ ZIP FOLDER (solution 1):               $(nadj_archive1).zip | downloading")
    run(`aws s3 cp s3://gda2020-ngca/adjustments/$(nadj_archive1).zip $(main_dir)/input/`)  
end
if isdir("$(nadj_dir1)")
    @info("$(now()) - NADJ FOLDER (solution 1):                   $(nadj_archive1)     | already exists")
else
    @info("$(now()) - NADJ FOLDER (solution 1):                   $(nadj_archive1)     | extracting")
    run(`unzip $(nadj_dir1).zip -d $(main_dir)/input/`)
end
# NADJ solution (2)
if isfile("$(nadj_dir2).zip")
    @info("$(now()) - NADJ ZIP FOLDER (solution 2):               $(nadj_archive2).zip | already exists")
else
    @info("$(now()) - NADJ ZIP FOLDER (solution 2):               $(nadj_archive2).zip | downloading")
    run(`aws s3 cp s3://gda2020-ngca/adjustments/$(nadj_archive2).zip $(main_dir)/input/`)  
end
if isdir("$(nadj_dir2)")
    @info("$(now()) - NADJ FOLDER (solution 2):                   $(nadj_archive2)     | already exists")
else
    @info("$(now()) - NADJ FOLDER (solution 2):                   $(nadj_archive2)     | extracting")
    run(`unzip $(nadj_dir2).zip -d $(main_dir)/input/`)
end


## ===
# ORGANSIE FILES
# - NADJ (solution 1)
# - NADJ (solution 2)
# - RVS

@info("$(now()) - READ FILES AND FORM DATAFRAMES:             $nadj_archive1     | $nadj_archive2 | $rvs_file")

# NADJ dataframe (solution 1)
# - Need to replace lon, lat because DynAdjust lon/lat is incorrect
df1 = dynxyz2dataframe(xyz_fp1)

# Replace lon/lat
tb_llh = columntable(xyz2llh.(df1.x, df1.y, df1.z))
lon = tb_llh[1]
lat = tb_llh[2]
df1.lon = lon
df1.lat = lat

# NADJ dataframe (solution 2)
# - Need to replace lon, lat because DynAdjust lon/lat is incorrect
df2 = dynxyz2dataframe(xyz_fp2)

# Replace lon/lat
tb_llh = columntable(xyz2llh.(df2.x, df2.y, df2.z))
lon = tb_llh[1]
lat = tb_llh[2]
df2.lon = lon
df2.lat = lat

# RVS dataframe
df_rvs = CSV.read(rvs_fp, DataFrame)
df_rvs = rename(df_rvs, :Column1 => :site)

## ===
# COORDINATE COMPARISON TO OTHER SOLUTION
# - NADJ(1) and NADJ(2)
@info("$(now()) - COORDINATE COMPARISON (NADJ):               $nadj_archive1     | $nadj_archive2")
df_nadjDelta = coordinate_delta_DYNADJUST(df1, df2)

## ===
# UNCERTAINTY COMPARISON TO OTHER SOLUTION
# - NADJ(1) and NADJ(2)
# - uncertainty_delta_SINEX() works for the NADJ dataframe too. 
@info("$(now()) - UNCERTAINTY COMPARISON (NADJ):              $nadj_archive1     | $nadj_archive2")
df_nadjUncDelta = uncertainty_delta_SINEX(df1, df2)

## ===
# COORDINATE COMPARISON TO RVS
@info("$(now()) - COORDINATE COMPARISON (RVS):                $nadj_archive1     | $rvs_file")
df_rvsDelta = coordinate_delta_to_RVS(df1, df_rvs)

## ===
# ADJUSTMENT SUMMARY BLOCK
# - from gda2020_YYYYMMDD.phased-stage.adj)

@info("$(now()) - READ ADJUSTMENT SUMMARY BLOCK:              gda2020_$(date1).phased-stage.adj")

# Read block
adjustment_summary_block = read_DynAdjustAdjustmentSummary_BLOCK("$(nadj_dir1)/gda2020_$(date1).phased-stage.adj")

# Number of measurements
i = findfirst(contains("Number of measurements"), adjustment_summary_block)
line = adjustment_summary_block[i]
num_measurements = split(line, " ", keepempty=false)[4]

# Number of outliers
i = findfirst(contains("Number of measurements"), adjustment_summary_block)
line = adjustment_summary_block[i]
num_outliers = split(line, " ", keepempty=false)[5]
num_outliers = chop(num_outliers, head=1, tail=0)

# Number of iterations
i = findlast(contains("ITERATION"), adjustment_summary_block)
line = adjustment_summary_block[i]
num_iterations = split(line, " ", keepempty=false)[2]

# Sigma zero
i = findfirst(contains("Sigma Zero"), adjustment_summary_block)
line = adjustment_summary_block[i]
sigma_zero = split(line, " ", keepempty=false)[4]

# Iteration time
i = findfirst(contains("Elapsed time"), adjustment_summary_block)
line = adjustment_summary_block[i]
iteration_time = split(line, " ", keepempty=false)[3]

# Convergence time
i = findfirst(contains("Total time"), adjustment_summary_block)
line = adjustment_summary_block[i]
convergence_time = split(line, " ", keepempty=false)[3]

## ===
# ADJUSTMENT COORDINATE CORRECTIONS 
# - from gda2020_YYYYMMDD.phased.cor

# Adjusted coordinates
df_coord_correction = dyncor2dataframe("$(nadj_dir1)/gda2020_$(date1).phased-stage.cor")

# Largest station shift
largest_station_shift = sort(df_coord_correction, order(:corr3D, rev=true))[1:1,:]

## ===
# LARGEST N-STATISTICS
# - from gda2020_YYYYMMDD.phased.adj
df_adjusted_measurements = dynadj2dataframe("$(nadj_dir1)/gda2020_$(date1).phased-stage.adj")

## ===
# MEASUREMENTS TO STATIONS
#df_m2s, df_m2s_warnings = dynm2s2dataframe("$(nadj_dir1)/gda2020_$(date1).m2s")

## ===
# COLLATE SUMMARY INFORMATION
# - Archive size
# - number of measurements
# - Additional sites
# - Removed sites
# - Input files

#  Archive size (MegaBytes)
zip_size1 = stat("$(nadj_dir1).zip").size/1000000
#zip_size2 = stat("$(nadj_dir2).zip").size/1000000

# Number of sites
num_sites1 = length(unique(df1.site))
num_sites2 = length(unique(df2.site))

# Non-common sites
additional_sites = setdiff(df1[:, "site"], df2[:, "site"])
removed_sites = setdiff(df2[:, "site"], df1[:, "site"])

# Input files
input_files_list_1 = sort(readdir("$nadj_dir1/inputFiles"))
input_files_list_2 = sort(readdir("$nadj_dir2/inputFiles"))
new_input_files = setdiff(input_files_list_1, input_files_list_2)
new_file_marker = []
for i in eachindex(input_files_list_1)
    if input_files_list_1[i] in new_input_files
        marker = "      (NEW INPUT)"
        push!(new_file_marker, marker)
    else
        marker = " "
        push!(new_file_marker, marker)
    end
end
df_input_files = DataFrame(
                            nadj1=input_files_list_1,
                            new=new_file_marker,
)

# Adjustment summary
@info("$(now()) - PRODUCE SUMMARY TEXT AND FIGURES.")

f = open("/$(output_dir)/GDA2020_NADJ_$(date1)_$(date2)_SUMMARY.txt", "w");

@printf(f, "\nHello all, \
            \n\nGDA2020(%s) is now available in the s3 bucket at gda2020-ngca/adjustments/.  \
            \n\nBelow are some details. \
            \n\nThanks, \
            \nGDA2020 Team\n\n", date1)

@printf(f, "\n--- SUMMARY ------------------------------------------------------------------------------------------\n\n")

@printf(f, "Solution name:                                     %s\n", nadj_archive1)
@printf(f, "Archive size:                                      %.2f MB\n", zip_size1)
@printf(f, "Number of sites:                                   %s\n", num_sites1)
@printf(f, "Number of measurements:                            %s\n", num_measurements)
@printf(f, "Number of outliers:                                %s\n", num_outliers)
@printf(f, "Number of iterations:                              %s\n", num_iterations)
@printf(f, "Iteration time (hours):                            ~%s\n", iteration_time)
@printf(f, "Convergence time (hours):                          %s\n", convergence_time)
@printf(f, "Sigma zero:                                        %s\n", sigma_zero)
@printf(f, "\n")
@printf(f, "Largest station shift (m):      %s\n", largest_station_shift.site[1])
@printf(f, "\n")
@printf(f, "                                %+7.5f (ΔE)\n" ,largest_station_shift.corrE[1]/1000)
@printf(f, "                                %+7.5f (ΔN)\n" ,largest_station_shift.corrN[1]/1000)
@printf(f, "                                %+7.5f (ΔU)\n" ,largest_station_shift.corrU[1]/1000)

@printf(f, "\n\n--- INPUT FILES --------------------------------------------------------------------------------------\n")
@printf(f, "* In comparison to %s\n\n", nadj_archive2)

pretty_table(f, 
             df_input_files, 
             tf=tf_borderless, 
             show_header=false,
             alignment=[:l, :r],
)


@printf(f, "\n\n--- ADJUSTMENT SUMMARY -------------------------------------------------------------------------------\n")
@printf(f, "* From DynAdjust solution\n\n")

for i in eachindex(adjustment_summary_block)
    @printf(f, "%s\n", adjustment_summary_block[i])
end

@printf(f, "\n\nTen largest coordinate corrections:\n")
@printf(f, "* Units: mm\n\n")

pretty_table(f, sort(df_coord_correction, order(:corr3D, rev=true))[1:10,:], show_subheader=false, formatters=ft_printf("%.2f"), tf=tf_borderless)

@printf(f, "\n--- TEN LARGEST N-STATS (NADJ) -----------------------------------------------------------------------\n")
@printf(f, "* From DynAdjust solution\n\n")

pretty_table(f, 
                 sort(select(df_adjusted_measurements, :type, :station1, :station2, :station3, :nstat), 
                 order(:nstat))[1:10,:], 
                 header=(["Type", "Station1", "Station2", "Station3", "N-STAT"]),  
                 formatters=ft_printf("%.2f", 5),
                 tf=tf_borderless,
                 alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n--- TEN LARGEST COORDINATE CHANGES (NADJ) ------------------------------------------------------------\n")
@printf(f, "* In comparison to %s (%s - %s)\n\n", nadj_archive2, date1, date2)

pretty_table(f, 
                 sort(select(df_nadjDelta, :site, :dE, :dN, :dU, :distance3D), 
                 order(:distance3D, rev=true))[1:10,:], 
                 header=(["Site", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                 ["           ", "(mm)", "(mm)", "(mm)", "(mm)"]), 
                 formatters=ft_printf("%.2f", [2,3,4,5]),
                 tf=tf_borderless,
                 alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n--- TEN LARGEST UNCERTAINTY CHANGES (NADJ) -----------------------------------------------------------\n")
@printf(f, "* In comparison to %s (%s - %s)\n\n", nadj_archive2, date1, date2)

pretty_table(f, 
                 sort(select(df_nadjUncDelta, :site, :dE, :dN, :dU, :distance3D), 
                 order(:distance3D, rev=true))[1:10,:], 
                 header=(["Site", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                 ["           ", "(mm)", "(mm)", "(mm)", "(mm)"]), 
                 formatters=ft_printf("%.2f", [2,3,4,5]),
                 tf=tf_borderless,
                 alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n--- TEN LARGEST COORDINATE CHANGES (RVS) -------------------------------------------------------------\n")
@printf(f, "* In comparison to 109 RVS stations (NADJ - RVS)\n\n")

pretty_table(f, 
                 sort(select(df_rvsDelta, :site, :dE, :dN, :dU, :distance3D), 
                 order(:distance3D, rev=true))[1:10,:], 
                 header=(["Site", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                 ["           ", "(mm)", "(mm)", "(mm)", "(mm)"]), 
                 formatters=ft_printf("%.2f", [2,3,4,5]),
                 tf=tf_borderless,
                 alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n--- STATION CHANGES ----------------------------------------------------------------------------------\n")
@printf(f, "* In comparison to %s\n\n", nadj_archive2)

@printf(f, "Sites added:                                       %s\n", length(additional_sites))
@printf(f, "Sites removed:                                     %s\n", length(removed_sites))
@printf(f, "\n")

# List station changes
if length(additional_sites) > 0
    @printf(f, "Sites added:\n")
    pretty_table(
                 f, 
                 additional_sites, 
                 tf=tf_borderless, 
                 show_header=false,
                 alignment=[:l],
    )
    @printf(f, "\n")
end
if length(removed_sites) > 0
    @printf(f, "Sites removed:\n")
    pretty_table(
                 f, 
                 removed_sites, 
                 tf=tf_borderless, 
                 show_header=false,
                 alignment=[:l],
    )
end

@printf(f, "\n\n------------------------------------------------------------------------------------------------------\n\n")

close(f)

## ===
# PLOT AND SAVE FIGURES

# NADJ-NADJ
fig_nadj_pdf = plot_delta_PDF(df_nadjDelta, "GDA2020 coordinate change between adjustment solutions ($date1 - $date2)")
fig_nadjUnc_pdf = plot_delta_PDF(df_nadjUncDelta, "GDA2020 coordinate change between adjustment solutions ($date1 - date2)")

savefig(fig_nadj_pdf, "$output_dir/coordinate_comparison_NADJ-$(date1)_NADJ-$(date2)_PDF.png")
savefig(fig_nadjUnc_pdf, "$output_dir/uncertainty_comparison_NADJ-$(date1)_NADJ-$(date2)_PDF.png")

# NADJ-RVS
fig_rvs_pdf = plot_delta_PDF(df_rvsDelta, "GDA2020 coordinate change between adjustment and RVS (RVS - $date1)")
fig_rvs_bar = plot_delta_BAR(df_rvsDelta, "GDA2020 coordinate change between adjustment and RVS (RVS - $date1)")

savefig(fig_rvs_pdf, "$output_dir/coordinate_comparison_NADJ-$(date1)_RVS_PDF.png")
savefig(fig_rvs_bar, "$output_dir/coordinate_comparison_NADJ-$(date1)_RVS_BAR.png")

@info("$(now()) - === NADJ SUMMARY COMPLETE ===")