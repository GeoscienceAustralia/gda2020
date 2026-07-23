import pandas as pd
import re

# Complete three checks
# 1. Check for duplicate solutions and that all solutions are present (ie if solution 1 and 3 are present, solution 2 should also be present)
# 2. Check that all disconts in the sinex file are in the disconts file
# 3. Check epochs format

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

def check_solutions(con_estimates_check, noncon_estimates_check):

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

    # Check for duplicates, remove 21NA and KALG as they are duped on purpose
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

    del df_noncon_estimates_check
    del df_con_estimates
    del df_all_estimates_check

def check_disconts(new_dis_seperate, new_aus_sol_estimate):
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

    del df_new_dis_seperate
    del df_new_aus_sol_estimates
    del solution_counts
    del sites_to_check
    del disc_pairs
    del site_pairs
    del missing

def check_epochs(new_aus_sol_estimates):

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

    del df_new_aus_sol_estimates