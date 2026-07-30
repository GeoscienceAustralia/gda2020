## ===
# DESCRIPTION
# - This script combines all the jurisdiction coordinate differences into
#   one total NGCA coordinate difference summary.
# - Useful for comparing whole NGCA processing tests. 
#
# Input:
#                       - Multiple Text files ([JUR]_NGCA_COORDINATE_DIFFERENCES_[DATE1]_[DATE2].txt).
#                       - Multiple DAT files ([JUR]__NGCA_YYYYMMDD_CLUSTERS.dat).
#
# Output:
#                       - Figure: coordinate_comparison_NGCA_[DATE1]_[DATE2]_PDF.png.
#                       - List: coordinate_comparison_NGCA_[DATE1]_[DATE2].txt.
#                       - List: clusters_NGCA_[DATE1].txt.

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
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/dynadjust_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/qaqc_operations.jl")
include("/home/ec2-user/REPOSITORY/Geodesy-Tools/geodetic_operations.jl")

# User input
dp = "/home/ec2-user/PROJECTS/GDA2020-QAQC/NGCA/PROD/output"
date1 = "20240715"
date2 = "20240715"

# ===
# COORDINATE DELTA DATAFRAME

# File list
file_list = filter(endswith("NGCA_COORDINATE_DIFFERENCES_$(date1)_$(date2).txt"), readdir("$dp/"))

# Read files
for i in eachindex(file_list)

    # Extract jurisdiction
    JUR = split(file_list[i], "_")[1]

    if i == 1
        global df_delta = CSV.read("$dp/$(file_list[i])", DataFrame, delim=' ', ignorerepeated=true, header=1, skipto=4)
        df_delta = insertcols(df_delta, 5, :jur=>fill(JUR, length(df_delta.soln)))
    else
        df_temp = CSV.read("$dp/$(file_list[i])", DataFrame, delim=' ', ignorerepeated=true, header=1, skipto=4)
        df_temp = insertcols(df_temp, 5, :jur=>fill(JUR, length(df_temp.soln)))
        df_delta = vcat(df_delta, df_temp)
    end
end

# Organise 
df_delta = sort(df_delta, order(:Δ3D, rev=true))
rename!(df_delta, :ΔE => :dE)
rename!(df_delta, :ΔN => :dN)
rename!(df_delta, :ΔU => :dU)
rename!(df_delta, :Δ3D => :d3D)

# ===
# CLUSTER DATAFRAME

# File list
file_list = filter(endswith("NGCA_$(date1)_CLUSTERS.dat"), readdir("$dp/"))

# Read files
for i in eachindex(file_list)

    # Extract jurisdiction
    JUR = split(file_list[i], "_")[1]

    if i == 1
        global df_clusters = CSV.read("$dp/$(file_list[i])", DataFrame, delim=' ', ignorerepeated=true, header=["cluster", "rinex", "antenna", "antHgt", "job"])
        df_clusters.jur = fill(JUR, length(df_clusters.cluster))
    else
        df_temp = CSV.read("$dp/$(file_list[i])", DataFrame, delim=' ', ignorerepeated=true, header=["cluster", "rinex", "antenna", "antHgt", "job"])
        df_temp.jur = fill(JUR, length(df_temp.cluster))
        df_clusters = vcat(df_clusters, df_temp)
    end

end

# ===
# SAVE

# Figure
fig = plot_delta_PDF(df_delta, "NGCA coordinate differences between September and July (AUSPOS-v3.0 and AUSPOS-v2.4)")
savefig(fig, "$dp/coordinate_comparison_NGCA_$(date1)_$(date2)_PDF.png")

# Text file
f = open("$dp/coordinate_comparison_NGCA_$(date1)_$(date2).txt", "w");
pretty_table(
                    f, 
                    sort(select(df_delta, :site, :soln, :cluster, :job, :jur, :dE, :dN, :dU, :d3D), 
                    order(:d3D, rev=true)),
                    header=(["site", "soln", "cluster", "job", "jur", "ΔE", "ΔN", "ΔU", "Δ3D"], 
                    [" ", " ", " ", " ", " ", "(mm)", "(mm)", "(mm)", "(mm)"]),
                    formatters=ft_printf("%.2f"), 
                    tf=tf_borderless,
                    alignment=[:l, :l, :l, :l, :l, :r, :r, :r, :r],
)
close(f)

# Text file
f = open("$dp/clusters_NGCA_$(date1).txt", "w");
pretty_table(
                    f, 
                    sort(select(df_clusters, :cluster, :rinex, :antenna, :antHgt, :job, :jur), 
                    order(:cluster)),
                    header=(["cluster", "rinex", "antenna", "antHgt", "job", "jur"]),
                    tf=tf_borderless,
                    alignment=[:l, :l, :l, :l, :l, :l],
)
close(f)
