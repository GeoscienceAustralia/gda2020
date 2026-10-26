#!/usr/bin/env python3
"""Translate conventional station IDs in DynaML station and measurement files.

By default, this preserves the directory conventions used by the original Perl
script:

* jurisdiction: parent directory of the current working directory
* translation tables: ~/transTables/<jurisdiction>TransTable*.csv
* APREF discontinuities: ~/apref/disconts????????.snx
* DynaML files: ~/ngca/<jurisdiction>/baselines/*_stn.xml

The special jurisdiction name ``surat`` is mapped to ``qld``.
"""

import argparse
import shutil
import tempfile
from pathlib import Path
from ngca_path import NGCA_DIR, NGCA_TRANS_TABLES_DIR, NGCA_APREF_DIR

parser = argparse.ArgumentParser()

parser.add_argument("-j", "--juris", required=True, help="Jurisdiction to process")

args = parser.parse_args()

# Directory
JURIS_DIR = NGCA_DIR / args.juris.lower()
SINEX_DIR = NGCA_DIR / args.juris.lower() / 'sinexFiles'
RINEXANTLS_DIR = NGCA_DIR / args.juris.lower() / 'rinexantls'
BASELINES_DIR = NGCA_DIR / args.juris.lower() / 'baselines'

STATION_TAGS = ("Name",)
MEASUREMENT_TAGS = ("First", "Second")
SINEX_CONTROL_MARKERS = ("*-----", "%ENDSNX", "-SOLUTION", "+SOLUTION")


def load_translation_tables(juris):
    """Load and copy matching translation tables.

    Lines whose first non-whitespace character is ``#`` are
    ignored. Conflicting duplicate source IDs are rejected rather than silently
    overwritten.
    """
    
    # Look for correct state trans table
    pattern = f"{juris}TransTable*.csv"
    tables = sorted(NGCA_TRANS_TABLES_DIR.glob(pattern))
    
    if not tables:
        raise FileNotFoundError(
            f"No translation tables match {NGCA_TRANS_TABLES_DIR / pattern}"
        )
    elif len(tables) > 1:
        raise RuntimeError(
            f"More then one trans table found matching {NGCA_TRANS_TABLES_DIR / pattern}"
        )
    else:
        table = tables[0]


    translations = {}

    # Iterate through each line of table and add to dictionary
    with table.open("r", encoding="utf-8-sig", newline="") as file:
        for line in file:
            line = line.split(",")
            
            # Skip if starts with #
            if line[0].lstrip().startswith("#"):
                continue

            # Raise error if not exactly 2 entries
            elif len(line) != 2:
                raise ValueError(
                    f"Malformed row in {table}, line {line}: "
                    "expected only two columns"
                )

            ngca_name = line[0].strip()
            state_name = line[1].strip()

            # Check if ngca_name has already been assigned
            if ngca_name in translations:
                if state_name != translations[ngca_name]:
                    raise ValueError(
                        f"Conflicting translations for {ngca_name!r}: "
                        f"{state_name} or {translations[ngca_name]}"
                    )
                continue

            translations[ngca_name] = state_name

    return translations, table


def load_apref_stations():
    """Read station IDs from APREF discontinuity SINEX files."""
    
    # Find APREF discont file
    files = sorted(NGCA_APREF_DIR.glob("disconts????????.snx"))
    stations = set()

    if not files:
        raise FileNotFoundError(
            f"No APREF disconts file found at {NGCA_APREF_DIR}"
        )
    elif len(files) > 1:
        raise RuntimeError(
            f"More  then one APREF disconts file found at {NGCA_APREF_DIR}"
        )
    else:
        discont_file = files[0]

    # Step through each line of discont file
    with discont_file.open("r", encoding="utf-8", errors="replace") as file:
        
        # Skip all lines until +SOLUTION/DISCONTINUITY flag
        sol_marker = "+SOLUTION/DISCONTINUITY"
        sol_tracker = False
        for line in file:
            if sol_marker in line:
                sol_tracker = True
                continue
            if sol_tracker:
                if any(marker in line for marker in SINEX_CONTROL_MARKERS):
                    continue
                fields = line.split()
                if fields:
                    stations.add(fields[0])

    return stations


def translate_xml_line(line, tags, translations, apref_stations,):
    """Translate the text value of the first matching XML element on a line."""
    
    
    for tag in tags:
        opening = f"<{tag}>"
        closing = f"</{tag}>"
        opening_index = line.find(opening)
        if opening_index == -1:
            continue

        value_start = opening_index + len(opening)
        value_end = line.find(closing, value_start)
        if value_end == -1:
            continue

        station_id = line[value_start:value_end]
        if not station_id or station_id[:4] in apref_stations:
            return line, None, None

        translated_name = translations.get(station_id)
        if translated_name is None or translated_name == station_id:
            return line, None, None

        translated_line = (
            line[:value_start] + translated_name + line[value_end:]
        )
        return translated_line, station_id, translated_name

    return line, None, None


def translate_xml_file(path, tags, translations, apref_stations,):
    """Translate one XML file atomically while preserving its line formatting."""
    replacement_count = 0
    replacements = set()

    # Open file in read mode so it doesnt get edited
    with path.open("r", encoding="utf-8", newline="") as source:
        
        # Create temp file
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            
            # Try to translate each line and add to temp file
            try:
                for line in source:
                    output_line, old_name, new_name = translate_xml_line(
                        line, tags, translations, apref_stations
                    )
                    temporary.write(output_line)
                    if old_name is not None and new_name is not None:
                        replacement_count += 1
                        replacements.add((old_name, new_name))
            except Exception:
                temporary_path.unlink(missing_ok=True)
                raise
    
    # replace file with temp file
    try:
        shutil.copymode(path, temporary_path)
        temporary_path.replace(path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    return replacement_count, replacements


def measurement_path_for(station_path):
    """Return the corresponding *_msr.xml path for a *_stn.xml file."""
    
    suffix = "_stn.xml"
    
    if not station_path.name.endswith(suffix):
        raise ValueError(f"Unexpected station filename: {station_path.name}")
    
    return station_path.with_name(station_path.name[: -len(suffix)] + "_msr.xml")


def main():

    jurisdiction = args.juris.lower()

    try:
        translations, table_file = load_translation_tables(jurisdiction)
        apref_stations = load_apref_stations()

        station_files = sorted(BASELINES_DIR.glob("*_stn.xml"))
        if not station_files:
            raise FileNotFoundError(
                f"No station files match {BASELINES_DIR / '*_stn.xml'}"
            )

        total_replacements = 0
        all_replacements = set()

        for station_file in station_files:
            measurement_file = measurement_path_for(station_file)
            if not measurement_file.is_file():
                raise FileNotFoundError(
                    f"Missing measurement file for {station_file.name}: "
                    f"{measurement_file}"
                )

            station_count, station_replacements = translate_xml_file(
                station_file, STATION_TAGS, translations, apref_stations
            )
            measurement_count, measurement_replacements = translate_xml_file(
                measurement_file, MEASUREMENT_TAGS, translations, apref_stations
            )

            file_count = station_count + measurement_count
            total_replacements += file_count
            all_replacements.update(station_replacements)
            all_replacements.update(measurement_replacements)


        print(f"Loaded {len(translations)} translations from {table_file}.")
        print(
            f"Translated {total_replacements} XML value(s) representing "
            f"{len(all_replacements)} unique station mapping(s)."
        )
        return 0

    except (OSError, ValueError) as error:
        print(f"Error: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
