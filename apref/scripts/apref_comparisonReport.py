from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from datetime import date
from geodepy import gnss as gnss
from geodepy import convert as convert
import sys
import glob
import pandas as pd
import math
import argparse

from gda2020.apref.scripts.apref_comparisonQAQC import check_solutions, check_disconts, check_epochs, check_pos
from gda2020.apref.scripts.apref_comparisonMap import make_map

"""
This script will create a transition report between different apref versions. 

Usage:
python3 apref_comparison_report.py --old {old_apref} --new {new_apref}

Input:
    old_apref: Old APREF versions (eg 20240713)
    new_apref: NEW APREF versions (eg 20260110)

    IMPORTANT
    All apref input files must be in a folder with the name entered above. The required files in the folder are as follows:
    
    Old APREF directory must contain:
      Full SINEX solution      : AUS0OPSSNX_*SOL.SNX
      AUS-only SINEX solution  : AUS0OPSSNX_*SOL.SNX.AUS
      APREF SINEX solution     : apref*.snx

    New APREF directory must contain:
      Full SINEX solution      : AUS0OPSSNX_*SOL.SNX
      AUS-only SINEX solution  : AUS0OPSSNX_*SOL.SNX.AUS
      APREF SINEX solution     : apref*.snx
      Non-constraint solution  : *.SNX.NONCON.AUS
      Discontinuities file     : disconts*.snx

Output: 
    APREF_Comparison_Report_{old_apref}_{new_apref}.pdf: PDF report showing APREF comparison
    apref_station_map.png: Map showing changes in APREF

"""

#------ Functions ------#

def add_location(df, location_list):
    """
    Add location description to a data frame with site names

    :param df: pandas.DataFrame with a "Site" column.
    :param location_list: list of tuples containing site information, where the first element is the site name and the fifth element is the location description.
    :return: pandas.DataFrame with an added "Location" column containing the state of each site.
    """

    site_location_lookup = {
        row[0].strip().upper(): row[4].strip()
        for row in location_list
    }

    # Add location description
    df["Location"] = (
        df["Site"]
        .map(site_location_lookup)
        .fillna("")  
    )

    # Normalise location
    df["Location"] = (
        df["Location"]
        .astype(str)
        .str.strip()
        .str.split(",")
        .str[0]
        .str.strip()
        .str.upper()
    )

    return df

def make_multi_column_table(df, columns, sort_by=None, n_blocks=3, repeated_headers=True):
    """
    Convert a long DataFrame into a wide multi-column table.

    Example with columns=["site", "location"]:
        Site | Location | Site | Location | Site | Location

    :param df: pandas.DataFrame to be made into columns
    :param columns: list of columns to include in each repeated block.
    :param sort_by: list of columns to sort by.
    :param n_blocks: Number of horizontal blocks across the page.
    :param repeated_headers: True, final headers repeat as Site, Location, Solution. False, headers become Site 1, Location 1, etc.
    :return: pandas.DataFrame with multiple columns.
    """

    df = df.copy()

    # Keep only required columns
    df = df[columns].copy()

    # Clean string columns
    for col in columns:
        df[col] = df[col].astype(str).str.strip()

    # Optional sort
    if sort_by is not None:
        df = df.sort_values(by=sort_by)

    rows_per_block = math.ceil(len(df) / n_blocks)

    blocks = []

    # Create blocks
    for i in range(n_blocks):
        block = df.iloc[
            i * rows_per_block:(i + 1) * rows_per_block
        ].reset_index(drop=True)

        if repeated_headers:
            # Headers repeat as Site, Location, Solution
            block.columns = columns
        else:
            # Headers become Site 1, Location 1, Solution 1, etc.
            block.columns = [f"{col} {i + 1}" for col in columns]

        blocks.append(block)

    wide_df = pd.concat(blocks, axis=1)

    return wide_df.fillna("")

# Remove trailing zeros from epoch
def format_epoch(epoch):
    """
    Remove trailing zeros from epoch

    :param epoch: Epoch string in YY:DOY:00000 format
    :return: Epoch string in YY:DOY format
    """
    return str(epoch).replace(":00000", "")

def get_unique_sites(sites_list):
    """
    Get unique sites from a list of sites.

    :param sites_list: List of tuples containing site information, where the first element is the site name.
    :return: List of unique sites.
    """
    unique_sites_names = []
    unique_sites = []
    for site in sites_list:
        if site[0] not in unique_sites_names:
            unique_sites_names.append(site[0])
            unique_sites.append(site)

    return unique_sites

def diff_sign_adder(diff):
    """
    Add a sign to the difference and return positive or negative for colouring.

    :param diff: The difference value.
    :return: A tuple containing the signed difference and its colour category.
    """
    if diff > 0:
        return "+" + str(diff), "positive"
    elif diff < 0:
        return "-" + str(abs(diff)), "negative"
    else:
        return str(diff), ""


#------ Setup ------#

env = Environment(loader=FileSystemLoader("templates"))
template = env.get_template("apref_report.html")

parser = argparse.ArgumentParser(
    description="Create report comparing APREF versions.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  script.py
  script.py -o 20240713 -n 20260110

Directory requirements
----------------------

Old APREF directory must contain:
  Full SINEX solution      : AUS0OPSSNX_*SOL.SNX
  AUS-only SINEX solution  : AUS0OPSSNX_*SOL.SNX.AUS
  APREF SINEX solution     : apref*.snx
  Non-constraint solution  : *.SNX.NONCON.AUS

New APREF directory must contain:
  Full SINEX solution      : AUS0OPSSNX_*SOL.SNX
  AUS-only SINEX solution  : AUS0OPSSNX_*SOL.SNX.AUS
  APREF SINEX solution     : apref*.snx
  Non-constraint solution  : *.SNX.NONCON.AUS
  Discontinuities file     : disconts*.snx

"""
)

parser.add_argument(
    "-o", "--old",
    required=False,
    help="Old APREF version (e.g. 20240713)"
)

parser.add_argument(
    "-n", "--new",
    required=False,
    help="New APREF version (e.g. 20260110)"
)

args = parser.parse_args()

print("Creating APREF Comparison Report...")

# Determine the dates of the old and the new solution
if args.old and args.new:
    old = args.old
    new = args.new
    if new < old:
        print('The date of the old solution is more recent than the date '
              + 'of the new solution')
        sys.exit()
else:
    dirs = glob.glob('20??????')
    dirs.sort()
    new = dirs[-1]
    old = dirs[-2]

old_apref_name = "APREF " + old
new_apref_name = "APREF " + new

print("  Old solution: " + old)
print("  New solution: " + new)

# Find all files needed
old_full_sol = None
for f in glob.glob(old + '/AUS0OPSSNX_*SOL.SNX'):
#for f in glob.glob(old + '/XVSOLFIN*.SNX'):   # Needed for older APREF naming
    old_full_sol = f
if old_full_sol is None:
    raise TypeError("Can't find old solution sinex file")

old_aus_sol = None
for f in glob.glob(old + '/AUS0OPSSNX_*_SOL.SNX.AUS'):
#for f in glob.glob(old + '/XVSOLFIN*.SNX.AUS'):   # Needed for older APREF naming
    old_aus_sol = f
if old_aus_sol is None:
    raise TypeError("Can't find old GDA2020 solution sinex file")

old_con_sol = None
for f in glob.glob(old + '/apref*.snx'):
    old_con_sol = f
if old_con_sol is None:
    raise TypeError("Can't find old constraining solution sinex file")

old_noncon_sol = None
for f in glob.glob(old + '/*.SNX.NONCON.AUS'):
    old_noncon_sol = f
if old_noncon_sol is None:
    raise TypeError("Can't find new non-constraining solution sinex file")

new_full_sol = None
for f in glob.glob(new + '/AUS0OPSSNX_*SOL.SNX'):
    new_full_sol = f
if new_full_sol is None:
    raise TypeError("Can't find new solution sinex file")

new_aus_sol = None
for f in glob.glob(new + '/AUS0OPSSNX_*_SOL.SNX.AUS'):
    new_aus_sol = f
if new_aus_sol is None:
    raise TypeError("Can't find new GDA2020 solution sinex file")

new_con_sol = None
for f in glob.glob(new + '/apref*.snx'):
    new_con_sol = f
if new_con_sol is None:
    raise TypeError("Can't find new constraining solution sinex file")

new_noncon_sol = None
for f in glob.glob(new + '/*.SNX.NONCON.AUS'):
    new_noncon_sol = f
if new_noncon_sol is None:
    raise TypeError("Can't find new non-constraining solution sinex file")

new_dis = None
for f in glob.glob(new + '/disconts*.snx'):
    new_dis = f
if new_dis is None:
    raise TypeError("Cant find new disontinuity file")

# Read in required sinex files
print("1. Reading sinex files...")
old_full_sol_sites = gnss.read_sinex_sites(old_full_sol)
old_aus_sol_sites = gnss.read_sinex_sites(old_aus_sol)
old_con_sol_sites = gnss.read_sinex_sites(old_con_sol)
new_full_sol_sites = gnss.read_sinex_sites(new_full_sol)
new_aus_sol_sites = gnss.read_sinex_sites(new_aus_sol)
new_con_sol_sites = gnss.read_sinex_sites(new_con_sol)

#------ Integrity Checks ------#

# 3 checks to be completed to ensure the integrity of the sinex file.

# 1. Check for duplicate solutions and that all solutions are present (ie if solution 1 and 3 are present, solution 2 should also be present)
# 2. Check that all disconts in the sinex file are in the disconts file
# 3. Check epochs format
# 4. Check how much sites have moved

print("2. Performing checks on sinex files...")

# 1. Check solutions
new_con_sol_estimates = gnss.read_sinex_estimate(new_con_sol)
new_noncon_sol_estimates = gnss.read_sinex_estimate(new_noncon_sol)
check_solutions(new_con_sol_estimates, new_noncon_sol_estimates)

# 2. Check disconts

# Read in discontinuty file and solutions in gda2020 area
new_dis_seperate = gnss.read_disconts(new_dis)
new_aus_sol_estimates = gnss.read_solution_epochs(new_aus_sol)

check_disconts(new_dis_seperate, new_aus_sol_estimates)

del new_dis_seperate

# 3. Check epoch format

check_epochs(new_aus_sol_estimates)

del new_aus_sol_estimates

# 4. Check moved sites

old_con_sol_estimates = gnss.read_sinex_estimate(old_con_sol)
old_noncon_sol_estimates = gnss.read_sinex_estimate(old_noncon_sol)

csv_name = f"APREF_Movement_{old}_{new}.csv"

check_pos(old_con_sol_estimates, new_con_sol_estimates, old_noncon_sol_estimates, new_noncon_sol_estimates, csv_name)

del old_con_sol_estimates
del new_con_sol_estimates
del old_noncon_sol_estimates
del new_noncon_sol_estimates

#------ Overview ------#

# Calculate differences in number of sites
total_diff = len(new_full_sol_sites) - len(old_full_sol_sites)
gda2020_diff = len(new_aus_sol_sites) - len(old_aus_sol_sites)

# Find unique sites
new_con_sol_sites = get_unique_sites(new_con_sol_sites)
old_con_sol_sites = get_unique_sites(old_con_sol_sites)

constraining_diff = len(new_con_sol_sites) - len(old_con_sol_sites)

total_diff_sign, total_diff_class = diff_sign_adder(total_diff)
gda2020_diff_sign, gda2020_diff_class = diff_sign_adder(gda2020_diff)
constraining_diff_sign, constraining_diff_class = diff_sign_adder(constraining_diff)


#------ Added/removed Stations ------#

print("3. Finding added and removed stations...")

# Create list of added and removed stations in constraining solution
added_stations = []
removed_stations = []

# Create list of stations names
names_new_con = {r[0] for r in new_con_sol_sites}
names_old_con = {r[0] for r in old_con_sol_sites}

# Create list of added and removed stations
added_stations = [(r[0], r[4]) for r in new_con_sol_sites if r[0] not in names_old_con]
removed_stations = [(r[0], r[4]) for r in old_con_sol_sites if r[0] not in names_new_con]

# Create dataframes for added and removed stations
added_stations_df = pd.DataFrame(added_stations, columns=["Site", "Location"])
removed_stations_df = pd.DataFrame(removed_stations, columns=["Site", "Location"])

# Clean up station names and areas
added_stations_df["Site"] = added_stations_df["Site"].str.strip().str.upper()
added_stations_df["Location"] = (
    added_stations_df["Location"]
    .astype(str)
    .str.strip()
    .str.split(",")
    .str[0]          # get the state only
    .str.strip()
    .str.upper()
)

removed_stations_df["Site"] = removed_stations_df["Site"].str.strip().str.upper()
removed_stations_df["Location"] = (
    removed_stations_df["Location"]
    .astype(str)
    .str.strip()
    .str.split(",")
    .str[0]          # get the state only
    .str.strip()
    .str.upper()
)

# Check length
added_stations_count = len(added_stations_df)
removed_station_count = len(removed_stations_df)

if (added_stations_count - removed_station_count) != constraining_diff:
    raise ValueError("Calculated station change is not the same as number of stations found to have changed")

# Use later for map
new_sites = set(added_stations_df["Site"])
removed_sites = set(removed_stations_df["Site"])

# Make data frame have columns
added_stations_df = make_multi_column_table(added_stations_df, ["Site", "Location"], sort_by=["Location", "Site"])
removed_stations_df = make_multi_column_table(removed_stations_df, ["Site", "Location"], sort_by=["Location", "Site"])


#------ New disconts ------#

print("4. Finding new discontinuities...")

old_dis = gnss.read_solution_epochs(old_con_sol)
new_dis = gnss.read_solution_epochs(new_con_sol)

dis_cols = [
    "Site",
    "code1",
    "Solution",
    "code2",
    "start_epoch",
    "end_epoch",
    "mean_epoch"
]

# Make dis arr a df
old_dis_df = pd.DataFrame(old_dis, columns=dis_cols)
new_dis_df = pd.DataFrame(new_dis, columns=dis_cols)

# Clean strings
for df in [old_dis_df, new_dis_df]:
    df["Site"] = df["Site"].astype(str).str.strip().str.upper()
    df["Solution"] = df["Solution"].astype(str).str.strip()

# Compare using site and solution
old_solution_keys = set(zip(old_dis_df["Site"], old_dis_df["Solution"]))

# Find only new disconts
df_new_disconts = new_dis_df[
    ~new_dis_df[["Site", "Solution"]]
    .apply(tuple, axis=1)
    .isin(old_solution_keys)
].copy()

# Remove disconts for new stations
added_stations_names = {row[0] for row in added_stations}

df_new_disconts = df_new_disconts[
    ~df_new_disconts["Site"].isin(added_stations_names)
].copy()

#Find number of new disconts
new_discontinuities_count = len(df_new_disconts)

df_new_disconts = add_location(df_new_disconts, new_con_sol_sites)

# Aggregate solutions from the same site
df_new_disconts = (
    df_new_disconts
    .groupby(["Site", "Location"], as_index=False)
    .agg({
        "Solution": lambda x: ", ".join(sorted(set(x), key=lambda v: int(v)))
    })
)

# Sort by location then site name
df_new_disconts = df_new_disconts.sort_values(by=["Location", "Site"])

# Used for map later
new_disconts_sites = set(df_new_disconts["Site"])

df_new_disconts = make_multi_column_table(df_new_disconts, ["Site", "Location", "Solution"])


#------ Changed Disconts ------#

print("5. Finding changed discontinuities...")

# Split in three sections
# 1. Removed disconts
# 2. Start Epoch changed
# 3. End Epoch changed (but not new discont)


# 1. Removed disconts

new_solution_keys = set(zip(new_dis_df["Site"], new_dis_df["Solution"]))

df_removed_disconts = old_dis_df[
    ~old_dis_df[["Site", "Solution"]]
    .apply(tuple, axis=1)
    .isin(new_solution_keys)
].copy()

# Remove stations that were removed
removed_stations_names = {row[0] for row in removed_stations}

df_removed_disconts = df_removed_disconts[
    ~df_removed_disconts["Site"].isin(removed_stations_names)
].copy()

df_removed_disconts = add_location(df_removed_disconts, new_con_sol_sites)

# Agregate Solutions
df_removed_disconts = (
    df_removed_disconts
    .groupby(["Site", "Location"], as_index=False)
    .agg({
        "Solution": lambda x: ", ".join(sorted(set(x), key=lambda v: int(v)))
    })
)

df_removed_disconts["Change"] = ("Removed")

# 2. Start Epoch changed

# Merge df so all epochs on same row
df_merged = pd.merge(
    old_dis_df,
    new_dis_df,
    on=["Site","code1", "code2", "Solution"],
    suffixes=("_old", "_new")
)

# Find solutions with changed start epoch
df_start_changed = df_merged[
    df_merged["start_epoch_old"] != df_merged["start_epoch_new"]
].copy()

# Remove stations that have been added or removed
df_start_changed = df_start_changed[
    ~df_start_changed["Site"].isin(added_stations_names)
].copy()

df_start_changed = df_start_changed[
    ~df_start_changed["Site"].isin(removed_stations_names)
].copy()

df_start_changed = add_location(df_start_changed, new_con_sol_sites)

# Add column to show what changed
df_start_changed["Change"] = (
    "Start: " + df_start_changed["start_epoch_old"].apply(format_epoch) + " → " + df_start_changed["start_epoch_new"].apply(format_epoch)
)

# Remove columns not needed and sort
df_start_changed = df_start_changed.drop(columns=["code1", "code2", "start_epoch_old", "end_epoch_old", "mean_epoch_old", "start_epoch_new", "end_epoch_new", "mean_epoch_new"])

# 3. End Epoch changed (but not new discont)

# Find solutions with changed end epoch
df_end_changed = df_merged[
    df_merged["end_epoch_old"] != df_merged["end_epoch_new"]
].copy()

# Remove stations that have been added or removed
df_end_changed = df_end_changed[
    ~df_end_changed["Site"].isin(added_stations_names)
].copy()

df_end_changed = df_end_changed[
    ~df_end_changed["Site"].isin(removed_stations_names)
].copy()

# Remove stations that have added disconts
df_end_changed = df_end_changed[
    ~df_end_changed["Site"].isin(new_disconts_sites)
]

# Remove stations that have no disconts after (This is just the end date extending to available data)
df_all_sol = pd.concat(
    [
       old_dis_df[["Site", "Solution"]],
       new_dis_df[["Site", "Solution"]] 
    ],
    ignore_index=True
).drop_duplicates()

df_all_sol["Solution"] = pd.to_numeric(df_all_sol["Solution"])
df_end_changed["Solution"] = pd.to_numeric(df_end_changed["Solution"])

all_solutions = set(zip(df_all_sol["Site"], df_all_sol["Solution"]))

df_end_changed["has_next_solution"] = df_end_changed.apply(
    lambda row: (row["Site"], row["Solution"] + 1) in all_solutions,
    axis=1
)

df_end_changed = df_end_changed[df_end_changed["has_next_solution"]].copy()

# Add location description
df_end_changed = add_location(df_end_changed, new_con_sol_sites)

# Add column to show what changed
df_end_changed["Change"] = (
    "End: " + df_end_changed["end_epoch_old"].apply(format_epoch) + " → " + df_end_changed["end_epoch_new"].apply(format_epoch)
)

# Remove columns not needed
df_end_changed = df_end_changed.drop(columns=["code1", "code2", "start_epoch_old", "end_epoch_old", "mean_epoch_old", "start_epoch_new", "end_epoch_new", "mean_epoch_new", "has_next_solution"])

report_cols = ["Site", "Location", "Solution", "Change"]

# Combine all changed disconts
df_all_changed = pd.concat(
    [
        df_end_changed[report_cols],
        df_start_changed[report_cols],
        df_removed_disconts[report_cols],
    ],
    ignore_index=True
)

changed_discontinuities_count = len(df_all_changed)

df_all_changed = make_multi_column_table(df_all_changed, report_cols, ["Location", "Site"], 2)


#------ Moved ------#

print("6. Finding moved stations...")

# Read movement csv created before
df_movement = pd.read_csv(csv_name)

# Remove stations that were non_constraint
df_movement = df_movement[
    df_movement["Solution_type"] != "Non-constraint"
]

df_movement = add_location(df_movement, new_con_sol_sites)

# Find stations that have moved by more then 10mm
df_east_changed = df_movement[
    abs(df_movement["East_moved"]) >= 0.010
].copy()

df_north_changed = df_movement[
    abs(df_movement["North_moved"]) >= 0.010
].copy()

df_height_changed = df_movement[
    abs(df_movement["Height_moved"]) >= 0.010
].copy()

# Create column to show movement
df_east_changed["Change"] = (
    "East Diff: "
    + df_east_changed["East_moved"].round(4).astype(str)
)

df_north_changed["Change"] = (
    "North Diff: "
    + df_north_changed["North_moved"].round(4).astype(str)
)

df_height_changed["Change"] = (
    "Height Diff: "
    + df_height_changed["Height_moved"].round(4).astype(str)
)

df_pos_changed = pd.concat(
    [
        df_east_changed[report_cols],
        df_north_changed[report_cols],
        df_height_changed[report_cols],
    ],
    ignore_index=True,
)

# Find how many individual stations for count
df_pos_unique = df_pos_changed.drop_duplicates(subset=["Site"])

moved_sites = set(df_pos_unique["Site"])

moved_stations_count = len(df_pos_unique)

# Find largest movement for overview
largest_movement_site = df_movement.iloc[0]["Site"]
largest_movement = round(df_movement.iloc[0]["Movement_3D"], 3)

df_pos_changed = make_multi_column_table(
    df_pos_changed,
    report_cols,
    ["Location", "Site"],
    2,
)


#------ Map ------#

print("7. Creating map of stations...")

estimate_cols = [
    "Site",
    "Solution",
    "refEpoch",
    "staX",
    "staY",
    "staZ",
    "staX_sd",
    "staY_sd",
    "staZ_sd",
]

#remove columns no longer needed
drop_cols = estimate_cols.copy()
drop_cols.pop(0)
drop_cols.pop(0)

con_estimates = gnss.read_sinex_estimate(new_con_sol)
noncon_estimates = gnss.read_sinex_estimate(new_noncon_sol)

df_con_estimates = pd.DataFrame(con_estimates, columns=estimate_cols)

# Have to add velocity columns to read in data frame then remove columns
for col in ["velX", "velY", "velZ", "velX_sd", "velY_sd", "velZ_sd"]:
    estimate_cols.append(col)
df_noncon_estimates = pd.DataFrame(noncon_estimates, columns=estimate_cols)
df_noncon_estimates = df_noncon_estimates.drop(columns=["velX", "velY", "velZ", "velX_sd", "velY_sd", "velZ_sd"])

df_all_estimates = pd.concat([df_con_estimates, df_noncon_estimates], ignore_index=True)

# Only have one row per site
df_all_estimates = df_all_estimates.drop_duplicates(subset=["Site"])

# Convert all estimates from xyz to llh
llh = [
    convert.xyz2llh(x, y, z)
    for x, y, z in zip(
        df_all_estimates["staX"],
        df_all_estimates["staY"],
        df_all_estimates["staZ"]
    )
]

df_all_estimates[["lat", "lon", "h"]] = pd.DataFrame(
    llh,
    index=df_all_estimates.index
)

df_all_estimates = df_all_estimates.drop(columns=drop_cols)

def classify(site):
    """
    Classify a site into a category for mapping.
    
    :param site: The site name to classify.
    :return: A str representing the category of the site.
    """
    if site in new_sites:
        return "New Station"
    elif site in removed_sites:
        return "Removed Station"
    elif site in new_disconts_sites:
        return "New Discontinuity"
    elif site in moved_sites:
        return "Moved > 10 mm"
    else:
        return "Existing Station"

df_all_estimates["Category"] = df_all_estimates["Site"].apply(classify)

# Generate map
make_map(df_all_estimates)


#------ Make final report ------#

print("8. Generating PDF report...")

author_name = "Kyran Cook"

if not author_name:
    print("Author name not set, using default 'National Geodesy Team'")
    author_name = "National Geodesy Team"

html = template.render(
    old_apref=old_apref_name,
    new_apref=new_apref_name,
    author=author_name,
    section="National Geodesy",
    generated_date=date.today().strftime("%d %B %Y"),
    total_old=len(old_full_sol_sites),
    total_new=len(new_full_sol_sites),
    total_diff=total_diff_sign,
    total_diff_class=total_diff_class,
    gda2020_old=len(old_aus_sol_sites),
    gda2020_new=len(new_aus_sol_sites),
    gda2020_diff=gda2020_diff_sign,
    gda2020_diff_class=gda2020_diff_class,
    constraining_old=len(old_con_sol_sites),
    constraining_new=len(new_con_sol_sites),
    constraining_diff=constraining_diff_sign,
    constraining_diff_class=constraining_diff_class,
    moved_gt_10mm_count=moved_stations_count,
    largest_2d_moved=largest_movement,
    largest_movement_site=largest_movement_site,
    apref_map_path="apref_station_map.png",
    added_stations_table=added_stations_df.to_html(index=False),
    added_stations_count=added_stations_count,
    removed_stations_table=removed_stations_df.to_html(index=False),
    removed_stations_count=removed_station_count,
    new_discontinuities_table=df_new_disconts.to_html(index=False),
    new_discontinuities_count=new_discontinuities_count,
    changed_discontinuities_table=df_all_changed.to_html(index=False),
    changed_discontinuities_count=changed_discontinuities_count,
    moved_stations_table=df_pos_changed.to_html(index=False),
    moved_stations_count=moved_stations_count,
)

report_name = f"APREF_Comparison_Report_{old}_{new}.pdf"

HTML(string=html, base_url=".").write_pdf(report_name)

print("Done!")