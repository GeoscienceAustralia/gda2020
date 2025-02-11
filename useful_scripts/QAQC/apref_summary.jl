## ===
# DESCRIPTION
# This script aims to summarise an APREF cumulative solution that
# has been transformed to GDA2020 and limited to the GDA2020 extent.
#
# Items summarised:     
#                       - SINEX file 
#                       - Discontinuity file 
#                       - Stations added / removed
#                       - Discontinuities added / removed
#                       - Coordinate comparison to GDA2020 RVS
#                       - Coordinate comparison to previous APREF solution
#                       - Coordinate comparison to IGS solution
#
# Input:                
#                       - aprefYYYYMMDD.snx (solution 1)
#                       - aprefYYYYMMDD.disconts (solution 1)
#                       - aprefYYYYMMDD.snx (solution 2)
#                       - aprefYYYYMMDD.disconts (solution 2)
#                       - RVS_GDA2020.txt
#                       - igsYYPWW_SNX.csv || IGS0OPSSNX_1994002_2024104_00U_CRD.csv               
#                       
# Output:           
#                       - Summary file (.txt)
#                       - Figures of coordinate deltas (Probability Density Functions)
#                       - Figures of coordinate deltas (Bar plots)
# ToDo:                 
#                       - Add Type-B Uncertainties to IGS SINEX file.
#                       - Uncertainty transformation.
#                       - List largest uncertainties.
#                       - Remove stations by area (i.e. restrict to GDA2020 extent). 

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
using StatsPlots

# Functions
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/sinex_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/geodetic_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/qaqc_operations.jl")

# User input
# - date of solution 1
# - date of solution 2
# - RVS file
# - IGS file (_CRD.SNX, _DSC.SNX are added to this)
# - working directory
# - Reference epoch of IGS file
# - IGS transformation condition (ITRF2020-->ITRF2014)
# - Save coordinate comparison lists to file (APREF, RVS, IGS)
date1 = "20240713"
date2 = "20220716"
rvs_file = "RVS_GDA2020.txt"
igs_file = "IGS0OPSSNX_1994002_2024216_00U"
main_dir = "/home/ec2-user/PROJECTS/GDA2020-QAQC/APREF/PROD"
igsRefEpoch = 2015.0
igsTransformation = true
saveCoodinateDeltaFiles = true

# Derived variables
sinex_file1 = "apref$(date1).snx"
sinex_file2 = "apref$(date2).snx"
disconts_file1 = "apref$(date1).disconts"
disconts_file2 = "apref$(date2).disconts"
snx_disconts_file1 = "disconts$(date1).snx"
snx_disconts_file2 = "disconts$(date2).snx"
input_dir = "$(main_dir)/input"
output_dir = "$(main_dir)/output"
igs_dir = "$(main_dir)/igs"
sinex_fp1 = "$(input_dir)/$(sinex_file1)"
sinex_fp2 = "$(input_dir)/$(sinex_file2)"
snx_disconts_fp1 = "$(input_dir)/$(snx_disconts_file1)"
snx_disconts_fp2 = "$(input_dir)/$(snx_disconts_file2)"
rvs_fp = "$(main_dir)/other/$(rvs_file)"
igs_sol_fp = "$(igs_dir)/$(igs_file)_CRD.SNX"
igs_dsc_fp = "$(igs_dir)/$(igs_file)_DSC.SNX"

# Logging
@info("$(now()) - === APREF SUMMARY BEGIN ===")
@info("$(now()) - INPUT SINEX FILE (1):                       $sinex_file1")  
@info("$(now()) - INPUT SINEX FILE (2):                       $sinex_file2")
@info("$(now()) - INPUT IGS:                                  $igs_file")
@info("$(now()) - INPUT RVS FILE:                             $rvs_file")
@info("$(now()) - INPUT SINEX DISCONTINUITY FILE (1):         $snx_disconts_file1")
@info("$(now()) - INPUT SINEX DISCONTINUITY FILE (2):         $snx_disconts_file2")

## ===
# FORM DATAFRAMES FROM INPUT FILES
# - APREF file 1 (SINEX)
# - APREF file 2 (SINEX)
# - IGS file (SINEX)
# - RVS file (txt)

# Dataframe
df1 = sinex2dataframe(sinex_fp1, false, false, "$snx_disconts_fp1")
df2 = sinex2dataframe(sinex_fp2, false, false, "$snx_disconts_fp2")
df_igs = sinex2dataframe(igs_sol_fp, false, false, igs_dsc_fp)
df_rvs = CSV.read(rvs_fp, DataFrame)

# Organise
df1 = sort(df1, order(:soln))
df2 = sort(df2, order(:soln))
df_igs = sort(df_igs, order(:soln))
df_rvs = rename(df_rvs, :Column1 => :site)

## ===
# COORDINATE TRANSFORMATION (IGS)
# - (1) ITRF2020@2015.0 --> ITRF2014@2015.0
# - (2) ITRF2014@RefEpoch --> ITRF2014@2020.0 (GDA2020)
# - Only make sense for sites on Australian plate

# Transformation (1)
if igsTransformation == true

    @info("$(now()) - APPLY TRANSFORMATION TO IGS FILE:           ITRF2020 --> ITRF2014 (@$(igsRefEpoch))")

    tb_igs_xyz_ITRF2014 = columntable(coordinate_transformation.(
                                                                    df_igs.x, 
                                                                    df_igs.y, 
                                                                    df_igs.z,
                                                                    "IERS", 
                                                                    -0.0014, 
                                                                    -0.0009,
                                                                    0.0014,
                                                                    -0.42,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    igsRefEpoch, 
                                                                    igsRefEpoch, 
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0, 
                                                                    0.0, 
                                                                    0.0)
    )

    # Replace columns
    df_igs.x = tb_igs_xyz_ITRF2014[1]
    df_igs.y = tb_igs_xyz_ITRF2014[2]
    df_igs.z = tb_igs_xyz_ITRF2014[3]

end

# Transformation (2)
@info("$(now()) - APPLY TRANSFORMATION TO IGS FILE:           ITRF2014@$(igsRefEpoch) --> ITRF2014@2020.0")

tb_igs_xyz_gda2020 = columntable(coordinate_transformation.(
                                                            df_igs.x, 
                                                            df_igs.y, 
                                                            df_igs.z,
                                                            "AGRS",
                                                            0.0, 
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            igsRefEpoch, 
                                                            2020.0, 
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            1.50379, 
                                                            1.18346, 
                                                            1.20716)
)

# Replace columns
df_igs.x = tb_igs_xyz_gda2020[1]
df_igs.y = tb_igs_xyz_gda2020[2]
df_igs.z = tb_igs_xyz_gda2020[3]

## === 
# UPDATE GEOGRAPHIC COORDINATES
tb_igs_llh_gda2020 = columntable(xyz2llh.(df_igs.x, df_igs.y, df_igs.z))
df_igs.lon = tb_igs_llh_gda2020[1]
df_igs.lat = tb_igs_llh_gda2020[2]
df_igs.ellHgt = tb_igs_llh_gda2020[3]

## ===
# APREF SINEX FILE SUMMARY
# - Number of stations
# - Number of solutions
# - Non-common stations
# - Non-common solutions

# Number of stations
num_sites1 = length(unique(df1.site))
num_sites2 = length(unique(df2.site))

# Number of solutions
num_solns1 = length(df1.soln)
num_solns2 = length(df2.soln)

# Non-common stations 
additional_sites = setdiff(df1[:, "site"], df2[:, "site"])
removed_sites = setdiff(df2[:, "site"], df1[:, "site"])

# Non-common solutions
# - Filter out 1st solution (_1900001)
# - Filter out single solution sites
additional_solns = setdiff(df1[:, "soln"], df2[:, "soln"])
additional_solns = filter(s -> !endswith(s, "_1900001"), additional_solns)
additional_solns = filter(s -> length(s)!=4, additional_solns)
removed_solns = setdiff(df2[:, "soln"], df1[:, "soln"])
removed_solns = filter(s -> !endswith(s, "_1900001"), removed_solns)
removed_solns = filter(s -> length(s)!=4, removed_solns)

## ===
# COORDINATE COMPARISON TO OTHER APREF SOLUTION
# - Dataframe
# - Figure

# Coordinate delta
df_aprefDelta = coordinate_delta_SINEX(df1, df2)

# Remove datatframes
df1_temp = nothing
df2_temp = nothing

## ===
# COORDINATE COMPARISON TO RVS
# - Dataframe
# - Figure
# - A subset of columns
# - Also add 3D distance

df_rvsDelta = coordinate_delta_to_RVS(df1, df_rvs)
df_rvsDeltaDisplay = select(df_rvsDelta, :soln, :site, :numSolns, :dE, :dN, :dU)
df_rvsDeltaDisplay = insertcols(df_rvsDeltaDisplay, :distance3D => sqrt.(df_rvsDeltaDisplay.dE.^2 + df_rvsDeltaDisplay.dN.^2 + df_rvsDeltaDisplay.dU.^2))

## ===
# COORDINATE COMPARISON TO IGS
# - Dataframe
# - Figure

df_igsDelta = coordinate_delta_SINEX(df1, df_igs)

## ===
# DISPLAY SUMMARY

@info("$(now()) - PRODUCE SUMMARY TEXT AND FIGURES.")

igs_file = split(igs_file, ".")[1]

f = open("/$(output_dir)/APREF_$(date1)_$(date2)_SUMMARY.txt", "w");

@printf(f, "\n-------------------------------------- SUMMARY ---------------------------------------\n\n")

@printf(f, "SINEX file:                                        %s\n", sinex_file1)
@printf(f, "SINEX file (source):                               N/A\n")
@printf(f, "Discontinuity file:                                %s\n", disconts_file1)
@printf(f, "Discontinuity file (source):                       %s\n", snx_disconts_file1)
@printf(f, "Number of sites:                                   %s\n", num_sites1)
@printf(f, "Number of site solutions:                          %s\n", num_solns1)
@printf(f, "\n")
@printf(f, "Statistics (coordinate difference to %s)\n", "apref$(date2)")
@printf(f, "* Units: mm\n")
@printf(f, "* N: %i \n\n", length(df_aprefDelta.dE))

df_tb = DataFrame()
df_tb = insertcols(df_tb, 1, :c => ["ΔE:", "ΔN:", "ΔU:"])
df_tb = insertcols(df_tb, 2, :min => [minimum(df_aprefDelta.dE), minimum(df_aprefDelta.dN), minimum(df_aprefDelta.dU)])
df_tb = insertcols(df_tb, 3, :max => [maximum(df_aprefDelta.dE), maximum(df_aprefDelta.dN), maximum(df_aprefDelta.dU)])
df_tb = insertcols(df_tb, 4, :mean => [mean(df_aprefDelta.dE), mean(df_aprefDelta.dN), mean(df_aprefDelta.dU)])
df_tb = insertcols(df_tb, 5, :std => [std(df_aprefDelta.dE), std(df_aprefDelta.dN), std(df_aprefDelta.dU)])
pretty_table(f, 
             df_tb,
             header=["         ", "Min     ", "Max     ", "Mean     ", "Std     "],
             formatters=ft_printf("%.2f    ", [2,3,4,5]),
             tf=tf_borderless,
             alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n\nStatistics (coordinate difference to RVS)\n")
@printf(f, "* Units: mm\n")
@printf(f, "* N: %i \n\n", length(df_rvsDelta.dE))

df_tb = DataFrame()
df_tb = insertcols(df_tb, 1, :c => ["ΔE:", "ΔN:", "ΔU:"])
df_tb = insertcols(df_tb, 2, :min => [minimum(df_rvsDelta.dE), minimum(df_rvsDelta.dN), minimum(df_rvsDelta.dU)])
df_tb = insertcols(df_tb, 3, :max => [maximum(df_rvsDelta.dE), maximum(df_rvsDelta.dN), maximum(df_rvsDelta.dU)])
df_tb = insertcols(df_tb, 4, :mean => [mean(df_rvsDelta.dE), mean(df_rvsDelta.dN), mean(df_rvsDelta.dU)])
df_tb = insertcols(df_tb, 5, :std => [std(df_rvsDelta.dE), std(df_rvsDelta.dN), std(df_rvsDelta.dU)])
pretty_table(f, 
             df_tb,
             header=["         ", "Min     ", "Max     ", "Mean     ", "Std     "],
             formatters=ft_printf("%.2f    ", [2,3,4,5]),
             tf=tf_borderless,
             alignment=[:l, :r, :r, :r, :r],
)

@printf(f, "\n\nStatistics (coordinate difference to IGS)\n")
@printf(f, "* Units: mm\n")
@printf(f, "* N: %i \n\n", length(df_igsDelta.dE))

df_tb = DataFrame()
df_tb = insertcols(df_tb, 1, :c => ["ΔE:", "ΔN:", "ΔU:"])
df_tb = insertcols(df_tb, 2, :min => [minimum(df_igsDelta.dE), minimum(df_igsDelta.dN), minimum(df_igsDelta.dU)])
df_tb = insertcols(df_tb, 3, :max => [maximum(df_igsDelta.dE), maximum(df_igsDelta.dN), maximum(df_igsDelta.dU)])
df_tb = insertcols(df_tb, 4, :mean => [mean(df_igsDelta.dE), mean(df_igsDelta.dN), mean(df_igsDelta.dU)])
df_tb = insertcols(df_tb, 5, :std => [std(df_igsDelta.dE), std(df_igsDelta.dN), std(df_igsDelta.dU)])
pretty_table(f, 
             df_tb,
             header=["         ", "Min     ", "Max     ", "Mean     ", "Std     "],
             formatters=ft_printf("%.2f    ", [2,3,4,5]),
             tf=tf_borderless,
             alignment=[:l, :r, :r, :r, :r],
)


@printf(f, "\n---------------------------------- STATION CHANGES -----------------------------------\n")
@printf(f, "* In comparison to %s \n  and %s\n\n", sinex_file2, snx_disconts_file2)

# Station, solution and discontinuity changes
@printf(f, "Stations added:                                    %s\n", length(additional_sites))
@printf(f, "Stations removed:                                  %s\n", length(removed_sites))
@printf(f, "Solutions added:                                   %s\n", length(additional_solns))
@printf(f, "Solutions removed:                                 %s\n", length(removed_solns))
@printf(f, "\n")
# List station changes
if length(additional_sites) > 0
    @printf(f, "Stations added:\n")
    pretty_table(f, additional_sites, tf=tf_borderless, show_header=false)
end
    if length(removed_sites) > 0
    @printf(f, "Stations removed:\n")
    pretty_table(f, removed_sites, tf=tf_borderless, show_header=false)
end
@printf(f, "\n")
# List solution
if length(additional_solns) > 0
    @printf(f, "Solutions added:\n")
    pretty_table(f, additional_solns, tf=tf_borderless, show_header=false)
end
if length(removed_solns) > 0
    @printf(f, "\nSolutions removed:\n")
    pretty_table(f, removed_solns, tf=tf_borderless, show_header=false)
end

@printf(f, "\n----------------------- TEN LARGEST COORDINATE CHANGES (APREF) -----------------------\n")
@printf(f, "* In comparison to %s\n\n", sinex_file2)

pretty_table(
              f, 
              sort(select(df_aprefDelta, :soln, :site, :numSolns, :dE, :dN, :dU, :distance3D), 
              order(:distance3D, rev=true))[1:10,:], 
              header=(["Station", "Station", "Count", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
              ["(solution ID)", "(code only)", "(Solutions)", "(mm)", "(mm)", "(mm)", "(mm)"]), 
              formatters=ft_printf("%.2f", [4,5,6,7]),
              tf=tf_borderless,
              alignment=[:l, :c, :c, :r, :r, :r, :r],
)

@printf(f, "\n------------------------ TEN LARGEST COORDINATE CHANGES (RVS) ------------------------\n")
@printf(f, "* In comparison to 109 RVS stations\n\n")

pretty_table(f, sort(df_rvsDeltaDisplay, order(:distance3D, rev=true))[1:10,:], 
                 header=(["Station", "Station", "Count", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
                 ["(solution ID)", "(name only)", "(Solutions)", "(mm)", "(mm)", "(mm)", "(mm)"]), 
                 formatters=ft_printf("%.2f", [4,5,6,7]),
                 tf=tf_borderless,
                 alignment=[:l, :c, :c, :r, :r, :r, :r],
)

@printf(f, "\n------------------------ TEN LARGEST COORDINATE CHANGES (IGS) ------------------------\n")
@printf(f, "* In comparison to IGS stations (%s_CRD.SNX)\n\n", igs_file)

pretty_table(
              f, 
              sort(select(df_igsDelta, :soln, :site, :numSolns, :dE, :dN, :dU, :distance3D), 
              order(:distance3D, rev=true))[1:10,:], 
              header=(["Station", "Station", "Count", "ΔE ", "ΔN ", "ΔU ", "Δ3D "], 
              ["(solution ID)", "(name only)", "(Solutions)", "(mm)", "(mm)", "(mm)", "(mm)"]), 
              formatters=ft_printf("%.2f", [4,5,6,7]),
              tf=tf_borderless,
              alignment=[:l, :c, :c, :r, :r, :r, :r],
)

@printf(f, "\n\n--------------------------------------------------------------------------------------\n\n")

close(f)

## ===
# PLOT FIGURES

# APREF-APREF
fig_apref_pdf = plot_delta_PDF(df_aprefDelta, "APREF coordinate change between solutions ($sinex_file1 - $sinex_file2)")
fig_apref_bar = plot_delta_BAR(df_aprefDelta, "APREF coordinate change between solutions ($sinex_file1 - $sinex_file2)")

savefig(fig_apref_pdf, "$output_dir/coordinate_comparison_APREF-$(date1)_APREF-$(date2)_PDF.png")
savefig(fig_apref_bar, "$output_dir/coordinate_comparison_APREF-$(date1)_APREF-$(date2)_BAR.png")

# APREF-RVS
fig_rvs_pdf = plot_delta_PDF(df_rvsDelta, "APREF coordinate change to RVS ($sinex_file1 - $rvs_file)")
fig_rvs_bar = plot_delta_BAR(df_rvsDelta, "APREF coordinate change to RVS ($sinex_file1 - $rvs_file)")

savefig(fig_rvs_pdf, "$output_dir/coordinate_comparison_APREF-$(date1)_RVS_PDF.png")
savefig(fig_rvs_bar, "$output_dir/coordinate_comparison_APREF-$(date1)_RVS_BAR.png")

# APREF-IGS|
fig_igs_pdf = plot_delta_PDF(df_igsDelta, "APREF coordinate change to IGS ($sinex_file1 - $(igs_file).SNX)")
fig_igs_bar = plot_delta_BAR(df_igsDelta, "APREF coordinate change to IGS ($sinex_file1 - $(igs_file).SNX)")

savefig(fig_igs_pdf, "$output_dir/coordinate_comparison_APREF-$(date1)_$(igs_file)_PDF.png")
savefig(fig_igs_bar, "$output_dir/coordinate_comparison_APREF-$(date1)_$(igs_file)_BAR.png")

## ===
# SAVE COODRINATE DELTA LIST TO FILES
# - Comparison to other APREF.
# - Comparison to RVS.
# - Comparison IGS cumulative solution.

if saveCoodinateDeltaFiles == true

    @info("$(now()) - PRODUCE COORDINATE DELTA TEXT FILES.")

    # Coordinate delta (APREF)
    f = open("$(output_dir)/coordinate_delta_apref$(date1)_APREF.txt", "w");

    @printf(f, "Coordinate Comparison:  %s - %s (largest --> smallest)\n\n", sinex_file1, sinex_file2)
    
    pretty_table(
                    f,
                    sort(select(df_aprefDelta, :site, :soln, :dE, :dN, :dU, :distance3D, :d2_std95_E, :d2_std95_N, :d2_std95_U, :delta_std95_E, :delta_std95_N, :delta_std95_U, :outlierE, :outlierN, :outlierU),
                    order(:distance3D, rev=true)), 
                    header=(["Site", "Soln", "ΔE", "ΔN", "ΔU", "Δ3D", "apref2_σ95_E", "apref2_σ95_N", "apref2_σ95_U", "ΔE_σ95", "ΔNσ95", "ΔU_σ95", "outlierE", "outlierN", "outlierU"],
                    [" ", " ", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(95%)", "(95%)", "(95%)"]),
                    formatters=ft_printf("%.2f", [3,4,5,6,7,8,9,10,11,12]),
                    tf=tf_borderless,
                    alignment=[:l, :l, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r],
    )
    close(f)

    # Corodinate delta (RVS)
    f = open("$(output_dir)/coordinate_delta_apref$(date1)_RVS.txt", "w");

    @printf(f, "Coordinate Comparison:  %s - %s (largest --> smallest)\n\n", sinex_file1, rvs_file)

    pretty_table(
                    f,
                    sort(select(df_rvsDelta, :site, :soln, :dE, :dN, :dU, :distance3D, :d2_std95_E, :d2_std95_N, :d2_std95_U, :delta_std95_E, :delta_std95_N, :delta_std95_U, :outlierE, :outlierN, :outlierU),
                    order(:distance3D, rev=true)),
                    header=(["Site", "Soln", "ΔE", "ΔN", "ΔU", "Δ3D", "rvs_σ95_E", "rvs_σ95_N", "rvs_σ95_U", "ΔE_σ95", "ΔN_σ95", "ΔU_σ95", "outlierE", "outlierN", "outlierU"],
                    [" ", " ", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(95%)", "(95%)", "(95%)"]),
                    formatters=ft_printf("%.2f", [3,4,5,6,7,8,9,10,11,12]),
                    tf=tf_borderless,
                    alignment=[:l, :l, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r],
    )
    close(f)

    # Coordinate delta (IGS)
    f = open("$(output_dir)/coordinate_delta_apref$(date1)_IGS.txt", "w");

    @printf(f, "Coordinate Comparison:  %s - %s.SNX (largest --> smallest)\n\n", sinex_file1, igs_file)

    pretty_table(
                    f,
                    sort(select(df_igsDelta, :site, :soln, :dE, :dN, :dU, :distance3D, :d2_std95_E, :d2_std95_N, :d2_std95_U, :delta_std95_E, :delta_std95_N, :delta_std95_U, :outlierE, :outlierN, :outlierU),
                    order(:distance3D, rev=true)),
                    header=(["Site", "Soln", "ΔE", "ΔN", "ΔU", "Δ3D", "igs_σ95_E", "igs_σ95_N", "igs_σ95_U", "ΔE_σ95", "ΔN_σ95", "ΔU_σ95", "outlierE", "outlierN", "outlierU"],
                    [" ", " ", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(mm)", "(95%)", "(95%)", "(95%)"]),
                    formatters=ft_printf("%.2f", [3,4,5,6,7,8,9,10,11,12]),
                    tf=tf_borderless,
                    alignment=[:l, :l, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r, :r],
    )
    close(f)

end

@info("$(now()) - === APREF SUMMARY COMPLETE ===")
