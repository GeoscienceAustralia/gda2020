#!/usr/bin/env python3
"""Run a minimally constrained adjustment for a jurisdiction's baseline clusters."""

import argparse
import re
import subprocess

from ngca_path import DYNADJUST_ENV, NGCA_APREF_DIR, NGCA_DIR

parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")

args = parser.parse_args()

# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
SINEX_DIR = NGCA_DIR / args.juris.lower() / 'sinexFiles'
RINEXANTLS_DIR = NGCA_DIR / args.juris.lower() / 'rinexantls'
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'


def get_apref_stations(apref_file):
    stations: set[str] = set()
    in_site_id = False

    with apref_file.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "-SITE/ID" in line:
                break
            if in_site_id:
                station = line[1:5]
                if station.strip():
                    stations.add(station)
            if "+SITE/ID" in line:
                in_site_id = True

    return stations


def run_command(command, log_file, cwd, env=DYNADJUST_ENV):
    """
    Log commands to be run to a log file and then run command

    :param command: List with command to be run
    :type command: List
    :param log_file: Path to log file
    :type log_file: pathlib.Path
    :param cwd: Location for command to be run
    :type cwd: pathlib.Path
    """
    log_file.write(subprocess.list2cmdline(command) + "\n")
    log_file.flush()
    subprocess.run(command, check=True, cwd=cwd, env=env)


def main():

    apref_files = sorted(NGCA_APREF_DIR.glob("apref*.snx"))
    discont_files = sorted(NGCA_APREF_DIR.glob("disconts*.snx"))

    if not apref_files:
        raise FileNotFoundError(f"No apref*.snx found in {JURIS_DIR}")
    if not discont_files:
        raise FileNotFoundError(f"No disconts*.snx found in {JURIS_DIR}")

    apref_file = apref_files[-1]
    discont_file = discont_files[-1]
    apref_stations = get_apref_stations(apref_file)

    # Find all stations files
    station_files = sorted(BASELINES_DIR.glob("*_all_stn.xml"))
    if not station_files:
        raise FileNotFoundError(f"No *_all_stn.xml files found in {BASELINES_DIR}")

    # Match station files with msr files
    for stn_file in station_files:
        msr_file = stn_file.with_name(stn_file.name.replace("_stn.xml", "_msr.xml"))
        if not msr_file.is_file():
            raise FileNotFoundError(
                f"Measurement file not found for {stn_file.name}: {msr_file.name}"
            )

        # Get network name
        network = stn_file.name.removesuffix("_stn.xml")
        stations = []
        constraints = []
        done = set()
        mac1_station = None

        # Find stations in file for command line entries
        with stn_file.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if "<Name>" not in line:
                    continue

                match = re.search(r"<Name>\s*(.*?)\s*</Name>", line)
                if not match:
                    continue

                stn = match.group(1)

                if re.fullmatch(r"MAC1_\d{7}", stn):
                    mac1_station = stn

                if stn in done or stn.startswith("MAC1"):
                    continue

                done.add(stn)
                stations.append(stn)

                if stn in apref_stations:
                    constraints.append(f"{stn},CCC")
                elif re.fullmatch(r"\w{4}_\d{7}", stn) and stn[:4] in apref_stations:
                    constraints.append(f"{stn},CCC")

        include_stations = ",".join(stations)
        constraint_string = ",".join(constraints)

        if not include_stations:
            raise RuntimeError(f"No stations found in {stn_file.name}")
        if not constraint_string:
            raise RuntimeError(f"No APREF constraints identified for {stn_file.name}")

        # Create log file
        log_path = BASELINES_DIR / network
        with log_path.open("w", encoding="utf-8") as log:

            # dnaimport
            
            cmd = [
                "dnaimport", "-n", network, "-r", "GDA2020",
                *apref_files,
                stn_file.name,
                msr_file.name,
                "--include-stns-assoc-msrs", include_stations,
                "--split-gnss-cluster-msrs",
            ]

            if mac1_station:
                cmd.extend(["--exclude-stns-assoc-msrs", mac1_station])

            cmd.extend([
                "--exclude-msr-types", "Y",
                "--discontinuity-file", discont_file.name,
            ])
            run_command(cmd, log, BASELINES_DIR)

            # dnareftran

            run_command(
                ["dnareftran", "-n", network, "-r", "GDA2020"],
                log,
                BASELINES_DIR,
            )

            # dnaadjust

            run_command(
                [
                    "dnaadjust", "-n", network, "--simult",
                    "--constraints", constraint_string,
                    "--output-adj-msr",
                    "--sort-adj-msr-field", "7",
                    "--stn-coord-types", "ENzPLHhXYZ",
                ],
                log,
                BASELINES_DIR,
            )


if __name__ == "__main__":
    main()
