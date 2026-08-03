#!/usr/bin/env python3

# This script is used to setup the file strucutre and files needed to complete each step of the GDA2020 adjustment
#
# It is split into 4 areas that can be set up individually or all together
# 1. APREF
# 2. NGCA
# 3. NADJ
# 4. QAQC

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import yaml

VALID_PROCESSES = ["APREF", "NGCA", "NADJ", "QAQC"]

def load_config(config_path):
    with open(config_path, "r") as file:
        return yaml.safe_load(file)

def parse_args():

    parser = argparse.ArgumentParser(
        description="Set up repository folders for GDA2020."
    )

    parser.add_argument(
        "-p",
        "--process",
        nargs="+",
        type=str.upper,
        choices=VALID_PROCESSES,
        help="One or more processes to set up. Choices: APREF, NGCA, NADJ, QAQC.",
    )

    parser.add_argument(
        "-a",
        "--all",
        action="store_true",
        help="Set up all processes.",
    )

    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root directory. Defaults to the current working directory.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be created without creating folders.",
    )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete selected process folders.",
    )

    return parser.parse_args()


def prompt_for_processes():
    """
    Interactively prompt the user to select which processes to set up.

    :return: List of selected processes.
    :rtype: list
    """
    # Print all options
    print("\nGDA2020 Repository Setup")
    print("------------------------------------")
    print("Select one or more processes to set up(or clean):")
    print()

    for index, process in enumerate(VALID_PROCESSES, start=1):
        print(f"{index}. {process}")

    print("5. ALL")
    print()

    # Prompt for user input with basic validation
    while True:
        user_input = input("Enter selection: ").strip()

        if not user_input:
            print("No selection made. Please try again.")
            continue

        selected = list_user_input(user_input)

        if selected:
            return selected

        print("Invalid selection. Please enter numbers, process names, or 'all'.")


def list_user_input(user_input):
    """
    Create a list of selected processes based on user input.

    :param user_input: User input string.
    :type user_input: str
    :return: List of selected processes.
    :rtype: list
    """
    # Incase user enters commas replace, then split string
    user_inputs = user_input.replace(",", " ").split()

    selected = []

    # Loop through each user input and find the corresponding process
    for input in user_inputs:
        value = input.strip().upper()

        # If all is selected, return all processes
        if value in {"ALL", "A", "5"}:
            return VALID_PROCESSES.copy()

        # If number is selected, find corresponding process
        if value.isdigit():
            index = int(value)

            if 1 <= index <= len(VALID_PROCESSES):
                process = VALID_PROCESSES[index - 1]
                selected.append(process)
            else:
                return []

        # If process name is selected, add to list
        elif value in VALID_PROCESSES:
            selected.append(value)

        # If any user input is wrong, return empty list to prompt user again
        else:
            return []

    # Remove duplicates while preserving order
    selected = list(dict.fromkeys(selected))

    return selected

def determine_processes(args):
    """
    Determine which processes to set up based on command-line arguments or user input.

    :param args: command-line arguments.
    :type args: argparse.Namespace
    :return: List of selected processes.
    :rtype: list
    """

    if args.all and args.process:
        raise SystemExit("Error: Use either --all or --process, not both.")

    if args.all:
        return VALID_PROCESSES.copy()

    if args.process:
        # remove duplicates
        process_list = list(dict.fromkeys(args.process))
        return process_list

    # If no processes are specified, prompt the user interactively
    return prompt_for_processes()

def setup_processes(repo_root, processes, dry_run=False):
    """
    Setup all specified processes.

    :param repo_root: Root directory of repo.
    :type repo_root: Path
    :param processes: List of processes to set up.
    :type processes: list
    :dry_run: If True, show what would be created without creating folders.
    :type dry_run: bool
    """
    repo_root = repo_root.resolve()

    print()
    print(f"Repository root: {repo_root}")
    print(f"Selected processes: {', '.join(processes)}")
    print()

    # Iterate through each process
    for process in processes:
        setup_single_process(
            repo_root=repo_root,
            process=process,
            dry_run=dry_run,
        )


def setup_single_process(repo_root, process, dry_run=False):
    """
    Setup a single process.

    :param repo_root: Root directory of repo.
    :type repo_root: Path
    :param process: Process to set up.
    :type process: str
    :param dry_run: If True, show what would be created without creating folders.
    :type dry_run: bool
    """

    print(f"Setting up {process}")

    # Setup folders for the process
    folders = CONFIG[process]["folders"]

    for folder in folders:
        folder_path = repo_root / folder

        if folder_path.exists():
            print(f"  Exists: {folder_path}")

        else:
            if dry_run:
                print(f"  Would create: {folder_path}")
            else:
                folder_path.mkdir(parents=True, exist_ok=True)
                print(f"  Created: {folder_path}")

    print()

def confirm_clean(processes):
    """
    Confirm with the user before cleaning specified processes.

    :param processes: List of processes to clean.
    :type processes: list
    """

    print("\nWARNING")
    print("The following processes will be premanently cleaned:")

    for process in processes:
        print(f"  {process}")

    print()
    response = input('Type "YES" to clean above processes:').strip()

    if response != "YES":

        raise SystemExit("Clean operation aborted by user.")

def clean_processes(repo_root, processes, dry_run=False):
    """
    Clean all specified processes.

    :param repo_root: Root directory of repo.
    :type repo_root: Path
    :param processes: List of processes to clean.
    :type processes: list
    :dry_run: If True, show what would be deleted without deleting folders.
    :type dry_run: bool
    """
    repo_root = repo_root.resolve()

    print()
    print(f"Repository root: {repo_root}")
    print(f"Selected processes to clean: {', '.join(processes)}")
    print()

    # Confirm with user before proceeding with deletion
    if not dry_run:
        confirm_clean(processes)

    # Iterate through each process
    for process in processes:
        for folder in CONFIG[process]["folders"]:
            folder_path = repo_root / folder

            if folder_path.exists():
                if dry_run:
                    print(f"  Would delete: {folder_path}")
                else:
                    # Remove the folder and its contents
                    shutil.rmtree(folder_path)
                    print(f"  Deleted: {folder_path}")
            else:
                print(f"  Does not exist: {folder_path}")
  

def main():
    args = parse_args()
    processes = determine_processes(args)
    config = load_config(args.repo_root / "config.yaml")

    if args.clean:
        clean_processes(
            repo_root=args.repo_root,
            processes=processes,
            dry_run=args.dry_run,
        )

    else:
        setup_processes(
            repo_root=args.repo_root,
            processes=processes,
            dry_run=args.dry_run,
        )

    print("Repository setup complete.")

# Load config file
CONFIG_FILE = Path(__file__).parent / "config.yaml"
CONFIG = load_config(CONFIG_FILE)
print(CONFIG)

if __name__ == "__main__":
    main()