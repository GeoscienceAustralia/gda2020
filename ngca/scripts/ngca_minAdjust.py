#!/usr/bin/env python3

"""
Run a minimally constrained adjustment for each baseline cluster.

Python conversion of ngca_minAdjust.pl.

20251111:
Include E N coordinates via --stn-coord-types, as XYZ files are retained.
"""

from pathlib import Path
import re
import shutil
import subprocess


HOME = Path.home()
WORK_DIR = Path.cwd()


def copy_required_files():
    """Copy DynaML schema, APREF SINEX, and discontinuity SINEX files."""

    # DynaML.xsd is only copied if it isn't already present.
    dynml = WORK_DIR / "DynaML.xsd"
    if not dynml.exists():
        shutil.copy2(HOME / "DynaML.xsd", dynml)

    # Remove existing APREF SINEX files.
    for path in WORK_DIR.glob("apref*.snx"):
        path.unlink()

    # Copy current APREF SINEX files.
    for src in (HOME / "apref").glob("apref*.snx"):
        shutil.copy2(src, WORK_DIR / src.name)

    # Remove existing discontinuity SINEX files.
    for path in WORK_DIR.glob("disconts*.snx"):
        path.unlink()

    # Copy current discontinuity SINEX files.
    for src in (HOME / "apref").glob("disconts*.snx"):
        shutil.copy2(src, WORK_DIR / src.name)


def get_apref_stations(apref_file):
    """Extract 4-character station names from the SINEX SITE/ID block."""

    stations = []
    in_site_id = False

    with apref_file.open() as fh:
        for line in fh:
            if "-SITE/ID" in line:
                break

            if in_site_id:
                stations.append(line[1:5])

            if "+SITE/ID" in line:
                in_site_id = True

    return stations


def run_command(command, log_file):
    """Log and execute a command."""

    command_text = " ".join(command)
    log_file.write(command_text + "\n")
    log_file.flush()

    subprocess.run(command, check=True)


def main():
    copy_required_files()

    apref_files = list(WORK_DIR.glob("apref*.snx"))
    discont_files = list(WORK_DIR.glob("disconts*.snx"))

    if not apref_files:
        raise FileNotFoundError("No apref*.snx file found")

    if not discont_files:
        raise FileNotFoundError("No disconts*.snx file found")

    # Matches the Perl script's effective behaviour of using one matching file.
    apref_file = apref_files[-1]
    discont_file = discont_files[-1]

    apref_stns = set(get_apref_stations(apref_file))

    # Loop over XML baseline station files.
    for stn_file in sorted(WORK_DIR.glob("*_all_stn.xml")):

        msr_file = Path(str(stn_file).replace("stn", "msr"))
        network = stn_file.name.replace("_stn.xml", "")

        stations = []
        constraints = []
        done = set()
        mac1_station = None

        with stn_file.open() as fh:
            for line in fh:

                if "<Name>" not in line:
                    continue

                match = re.search(r"<Name>\s*(.*?)\s*</Name>", line)
                if not match:
                    continue

                stn = match.group(1)

                # MAC1 discontinuity station.
                if re.fullmatch(r"MAC1_\d{7}", stn):
                    mac1_station = stn

                # Perl excludes all MAC1 names from the include list.
                if stn in done or stn.startswith("MAC1"):
                    continue

                done.add(stn)
                stations.append(stn)

                # Direct APREF station match.
                if stn in apref_stns:
                    constraints.append(f"{stn},CCC")
                    continue

                # Discontinuity station such as ABCD_2025001.
                if re.fullmatch(r"\w{4}_\d{7}", stn):
                    root = stn[:4]

                    if root in apref_stns:
                        constraints.append(f"{stn},CCC")

        include_stations = ",".join(stations)
        constraint_string = ",".join(constraints)

        # Command log, retaining the behaviour of the Perl script.
        with open(network, "w") as log:

            #
            # dnaimport
            #

            cmd = [
                "dnaimport",
                "-n", network,
                "-r", "GDA2020",
                *[str(p) for p in sorted(WORK_DIR.glob("apref*.snx"))],
                str(stn_file),
                str(msr_file),
                "--include-stns-assoc-msrs",
                include_stations,
                "--split-gnss-cluster-msrs",
            ]

            if mac1_station:
                cmd.extend([
                    "--exclude-stns-assoc-msrs",
                    mac1_station,
                ])

            cmd.extend([
                "--exclude-msr-types", "Y",
                "--discontinuity-file", str(discont_file),
            ])

            run_command(cmd, log)

            #
            # dnareftran
            #

            cmd = [
                "dnareftran",
                "-n", network,
                "-r", "GDA2020",
            ]

            run_command(cmd, log)

            #
            # dnaadjust
            #

            cmd = [
                "dnaadjust",
                "-n", network,
                "--simult",
                "--constraints", constraint_string,
                "--output-adj-msr",
                "--sort-adj-msr-field", "7",
                "--stn-coord-types", "ENzPLHhXYZ",
            ]

            run_command(cmd, log)


if __name__ == "__main__":
    main()
