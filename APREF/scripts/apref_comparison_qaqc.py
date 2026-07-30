import pandas as pd
import re
import numpy as np

from geodepy import convert as convert

# Complete four checks
# 1. Check for duplicate solutions and that all solutions are present (ie if solution 1 and 3 are present, solution 2 should also be present)
# 2. Check that all disconts in the sinex file are in the disconts file
# 3. Check epochs format
# 4. Check for sites that have moved significantly

ESTIMATE_COLS = [
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

def check_solutions(con_estimates_check, noncon_estimates_check):
    """
    Checks there are no duplicate solutions numbers and that the remaining are correct

    :param con_estimates_check: Array with new constraint site estimates using read_sinex_estimate
    :param noncon_estimates_check: Array with new non-constraint site estimates using read_sinex_estimate
    """

    estimate_cols = ESTIMATE_COLS.copy()

    vel_cols = ["velX", "velY", "velZ", "velX_sd", "velY_sd", "velZ_sd"]

    df_con_estimates = pd.DataFrame(con_estimates_check, columns=estimate_cols)

    # Have to add velocity columns to read in data frame then remove columns
    for col in vel_cols:
        estimate_cols.append(col)
    df_noncon_estimates_check = pd.DataFrame(noncon_estimates_check, columns=estimate_cols)
    df_noncon_estimates_check = df_noncon_estimates_check.drop(columns=vel_cols)

    for col in vel_cols:
        estimate_cols.remove(col)

    df_all_estimates_check = pd.concat([df_con_estimates, df_noncon_estimates_check], ignore_index=True)

    # Check for duplicates, remove ALIC and MOBS as they are duped on purpose
    duplicates = df_all_estimates_check[df_all_estimates_check.duplicated(subset=["Site", "Solution"], keep=False)].query("Site not in ['ALIC', 'MOBS']").sort_values(by=["Site", "Solution"])

    missing_solutions = {}

    for site, group in df_all_estimates_check.groupby("Site"):
        sols = sorted(group["Solution"].unique().astype(int))

        expected = set(range(min(sols), max(sols) + 1))
        missing = sorted(expected - set(sols))

        if missing:
            missing_solutions[site] = missing

    if len(duplicates) > 0:
        print("  2.1 Duplicates solutions found:")
        print(duplicates)
    else:
        print("  2.1 No duplicate solutions found")
    if len(missing_solutions) > 0:
        print("  2.2 Missing solutions found:")
        for site, solutions in missing_solutions.items():
            print(f"    {site}: {solutions}")
    else:
        print("  2.2 No missing solutions found")

def check_disconts(new_dis_seperate, new_aus_sol_estimate):
    """
    Checks for duplicate, incorrect or missing disconts

    :param new_dis_seperate: Array of discontinuties from the discontsYYYYMMDD.snx file using read_disconts
    :param new_aus_sol_estimates: Array with new site estimates for all sites within GDA2020 using read_sinex_estimate
    :return: Prints disconts that are incorrect
    """
    dis_cols = [
        "Site",
        "Code1",
        "Solution",
        "Code2",
        "start_epoch",
        "end_epoch",
        "Type",
    ]

    df_new_dis_seperate = pd.DataFrame(new_dis_seperate, columns=dis_cols)
    df_new_aus_sol_estimates = pd.DataFrame(new_aus_sol_estimate, columns=["Site","Point","Solution","obs","start","end","mean"])

    # Remove V type disconts
    df_new_dis_seperate = df_new_dis_seperate[df_new_dis_seperate["Type"] == "P"]

    # Make list of stations that have multiple solutions (as single solutions sites not in discont file)
    solution_counts = (df_new_aus_sol_estimates.groupby("Site")["Solution"].nunique())
    multi_solution_sites = solution_counts[solution_counts > 1].index

    # Create new df of just these sites
    sites_to_check = df_new_aus_sol_estimates[df_new_aus_sol_estimates["Site"].isin(multi_solution_sites)]

    # Create sets of both
    disc_pairs = set(zip(df_new_dis_seperate["Site"], df_new_dis_seperate["Solution"]))
    site_pairs = set(zip(sites_to_check["Site"], sites_to_check["Solution"]))

    # Find solutions not in discont file
    missing = site_pairs - disc_pairs

    if len(missing) > 0:
        print("  2.3 Found discontinuities not in discont file")
        for site, point in sorted(missing):
            print(f"    Site: {site} Point: {point}")
    else:
        print("  2.3 No discontinuities are missing from discont file")

def check_epochs(new_aus_sol_estimates):
    """
    Check epochs in solutions file are formatted correctly

    :param new_aus_sol_estimates: Array with new site estimates for all sites within GDA2020 using read_sinex_estimate
    :return: Prints epochs that are malformed
    """

    df_new_aus_sol_estimates = pd.DataFrame(new_aus_sol_estimates, columns=["Site","Point","Solution","obs","start","end","mean"])

    epoch_cols = ["start", "end", "mean"]

    epoch_pattern = re.compile(r"^\d{2}:\d{3}:\d{5}$")

    bad_epochs = []

    for _, row in df_new_aus_sol_estimates.iterrows():
        for col in epoch_cols:
            epoch = str(row[col]).strip()

            problem = None

            # Check format
            if not epoch_pattern.match(epoch):
                problem = "Invalid format"

            else:
                yy, doy, sec = epoch.split(":")
                doy = int(doy)
                sec = int(sec)

                if not (1 <= doy <= 366):
                    problem = f"Invalid DOY ({doy})"

                elif not (0 <= sec <= 86399):
                    problem = f"Invalid seconds ({sec})"

            if problem:
                bad_epochs.append({
                    "Site": row["Site"],
                    "Solution": row["Solution"],
                    "Column": col,
                    "Epoch": epoch,
                    "Problem": problem
                })

    df_bad_epochs = pd.DataFrame(bad_epochs)

    if len(df_bad_epochs) > 0:
        print("  2.4 Found incorrect epochs:")
        for _, row in df_bad_epochs.iterrows():
            print(f"    Site: {row['Site']}  Epoch: {row['Epoch']}  Problem: {row['Problem']}")
    else:
        print("  2.4 No incorrect epochs found")

def xyz2enu(row):
    """
    Convert XYZ corrdinates in dataframe to ENU

    :param row: pandas.Series with "staX", "staY", "staZ" columns.
    :return: pandas.Series with "east", "north", "h", "zone"
    """
    lat, lon, h = convert.xyz2llh(row["staX"], row["staY"], row["staZ"])
    _, zone, east, north, _, _ = convert.geo2grid(lat, lon)
    return pd.Series([east, north, h, zone])

def check_pos(
    old_con_estimate,
    new_con_estimate,
    old_noncon_estimate,
    new_noncon_estimate,
    output_csv,
):
    """
    Create a CSV showing station movement between two APREF versions
    for both constraint and non-constraint solutions.

    :param old_con_estimate: Array with old constraint site estimates using read_sinex_estimate
    :param new_con_estimate: Array with new constraint site estimates using read_sinex_estimate
    :param old_noncon_estimate: Array with old non-constraint site estimates using read_sinex_estimate
    :param new_noncon_estimate: Array with new non-constraint site estimates using read_sinex_estimate
    :return: Prints 10 largest movements and creates csv of all site movements
    """

    estimate_cols = ESTIMATE_COLS.copy()

    old_con_estimate_df = pd.DataFrame(old_con_estimate, columns=estimate_cols)
    new_con_estimate_df = pd.DataFrame(new_con_estimate, columns=estimate_cols)

    vel_cols = ["velX", "velY", "velZ", "velX_sd", "velY_sd", "velZ_sd"]

    # Have to add velocity columns for noncon files
    for col in vel_cols:
        estimate_cols.append(col)

    old_noncon_estimate_df = pd.DataFrame(old_noncon_estimate, columns=estimate_cols)
    old_noncon_estimate_df = old_noncon_estimate_df.drop(columns=vel_cols)

    new_noncon_estimate_df = pd.DataFrame(new_noncon_estimate, columns=estimate_cols)
    new_noncon_estimate_df = new_noncon_estimate_df.drop(columns=vel_cols)   

    for col in vel_cols:
        estimate_cols.remove(col)

    # Convert all coords to enu
    old_con_estimate_df[["East", "North", "Height", "Zone"]] = old_con_estimate_df.apply(xyz2enu, axis=1)
    new_con_estimate_df[["East", "North", "Height", "Zone"]] = new_con_estimate_df.apply(xyz2enu, axis=1)
    old_noncon_estimate_df[["East", "North", "Height", "Zone"]] = old_noncon_estimate_df.apply(xyz2enu, axis=1)
    new_noncon_estimate_df[["East", "North", "Height", "Zone"]] = new_noncon_estimate_df.apply(xyz2enu, axis=1)

    # Remove columns no longer needed
    drop_cols = estimate_cols.copy()
    drop_cols.pop(0)
    drop_cols.pop(0)

    old_con_estimate_df = old_con_estimate_df.drop(columns=drop_cols)
    new_con_estimate_df = new_con_estimate_df.drop(columns=drop_cols)

    old_noncon_estimate_df = old_noncon_estimate_df.drop(columns=drop_cols)
    new_noncon_estimate_df = new_noncon_estimate_df.drop(columns=drop_cols)

    # Merge constraint then add tag
    merged_con = pd.merge(old_con_estimate_df, new_con_estimate_df, on=["Site", "Solution"], suffixes=("_old", "_new"))

    merged_con["Solution_type"] = "Constraint"

    # Merge non-constraint then add flag
    merged_noncon = pd.merge(old_noncon_estimate_df, new_noncon_estimate_df, on=["Site", "Solution"], suffixes=("_old", "_new"))

    merged_noncon["Solution_type"] = "Non-constraint"

    merged_df = pd.concat([merged_con, merged_noncon], ignore_index=True)

    # Calculate amount moved
    merged_df["East_moved"] = round(merged_df["East_new"] - merged_df["East_old"], 4)
    merged_df["North_moved"] = round(merged_df["North_new"] - merged_df["North_old"], 4)
    merged_df["Height_moved"] = round(merged_df["Height_new"] - merged_df["Height_old"], 4)

    merged_df["Movement_3D"] = round(np.sqrt(
        merged_df["East_moved"] ** 2
        + merged_df["North_moved"] ** 2
        + merged_df["Height_moved"] ** 2
    ), 4)

    # Sort by largest movement
    merged_df = merged_df.sort_values("Movement_3D", ascending=False)

    # Keep only columns needed
    movement_df = merged_df[["Site", "Solution", "Solution_type", "East_moved", "North_moved", "Height_moved", "Movement_3D"]]

    print("  2.5 List of the 10 largest site movements: ")
    print(movement_df.head(10))

    # Write CSV
    merged_df.to_csv(output_csv, index=False)
