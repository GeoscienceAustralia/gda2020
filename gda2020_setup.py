#!/usr/bin/env python3

# This script is used to setup the file strucutre and files needed to complete each step of the GDA2020 adjustment
#
# It is split into 4 areas that can be set up individually or all together
# 1. APREF
# 2. NGCA
# 3. NADJ
# 4. QAQC

import argparse
from pathlib import Path
import shutil
import subprocess
import os

if shutil.which("pip") is None:
    subprocess.run(["sudo", "apt", "update"], check=True)
    subprocess.run(["sudo", "apt", "install", "-y", "python3-pip"], check=True)

REQUIRED_PACKAGES = {
    "requests": "requests",
    "PyYAML": "yaml",
}

# Install required packages if needed
for package, module in REQUIRED_PACKAGES.items():
    try:
        __import__(module)
    except ImportError:
        print("Installing requests")
        subprocess.run(["pip", "install", package], check=True)

import yaml
import requests

VALID_PROCESSES = ["APREF", "NGCA", "NADJ", "QAQC"]

JURIS = ["ACT", "NSW", "NT", "QLD", "SA", "TAS", "VIC", "WA"]

def load_config(config_path):
    # Load config file and validate it contains correct number of processes
    try:
        with open(config_path, "r") as file:
            config = yaml.safe_load(file)

            required_sections = set(VALID_PROCESSES + ["metadata"])

            # Check if config has exactly sections required
            if set(config.keys()) != required_sections:
                raise ValueError(f"Config file must contain exactly {sorted(required_sections)}")
            return config

    # Handle file not found and YAML parsing errors
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found: {config_path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML file: {e}")

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
    for item in user_inputs:
        value = item.strip().upper()

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
    print("Downloading required packages")
    print()

    # Downlad packages needed for rest of setup to run.
    for package in ["zip", "unzip"]:
        subprocess.run(["sudo", "apt", "install", "-y", package], check=True)
        print(f"  Downloaded {package}")

    if shutil.which("aws") is None:
        
        home_dir = os.path.expanduser("~")
        zip_file = os.path.join(home_dir, "awscliv2.zip")
        print(zip_file)
        
        subprocess.run(["curl", "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip", "-o", zip_file], check=True) # Using a link here is a bit dangerous but it is what aws recommends https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
        subprocess.run(["unzip", zip_file], cwd=home_dir, check=True)
        os.remove(zip_file)
        subprocess.run(["sudo", os.path.join(home_dir, "aws", "install")], check=True)
        print("  Downloaded AWS CLI")
    else:
        print("  AWS CLI already downloaded")

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

    print()
    print("Creating folders:")
    # Setup folders for the process
    folders = CONFIG[process].get("folders", [])

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
    print("Downloading files:")

    # Download files for process
    downloads = CONFIG[process].get("downloads", [])

    for download in downloads:

        # If the download is jurisdiction dependant, loop through each jurisdiction and download files for each one
        if download.get("jurisdiction_dependant"):
            for jurisdiction in JURIS:

                # If the jurisdiction_upper_case flag is not set or is False, convert the jurisdiction to lower case
                if not download.get("jurisdiction_upper_case", False):
                    jurisdiction = jurisdiction.lower()

                # If the include pattern is not specified, raise an error
                if not download.get("include"):
                    raise ValueError(f"  Download for {process} with jurisdiction_dependant=True must have an 'include' pattern specified.")

                # Download the files for the jurisdiction
                aws_download(uri=download["uri"], 
                             destination=download["destination"], 
                             include=f"{jurisdiction}{download['include']}", 
                             recursive=download.get("recursive", False), 
                             most_recent=download.get("most_recent", False),
                             preserve_folder=download.get("preserve_folder", False),
                             dry_run=dry_run)

        # Otherwise, download the files as specified in the config
        else:
            aws_download(**download, dry_run=dry_run)

        # Unzip the downloaded files if unzip flag is set to True
        if download.get("unzip", False):
            # Unzip the downloaded files if unzip flag is set to True
            destination_path = Path(download["destination"])
            if not destination_path.exists():
                raise ValueError(f"  Destination path does not exist: {destination_path}")

            for zip_file in destination_path.glob("*.zip"):
                if dry_run:
                    print(f"  Would unzip: {zip_file}")
                else:
                    print(f"  Unzipping {zip_file}")
                    shutil.unpack_archive(zip_file, destination_path)
                    print(f"  Unzipped: {zip_file}")
                    Path.unlink(zip_file)

    print()
    print("Downloading packages:")

    packages = CONFIG[process].get("packages", [])

    download_packages(packages, dry_run)

    print()

def find_latest_by_name(uri, destination, preserve_folder=False):
    """
    Find the most recent file or folder at an S3 URI path and return the updated URI and destination path.
    This is completed by looking at the name of file/folder so using YYYMMDD in naming is advised.

    :param uri: S3 URI to search
    :type uri: str
    :param destination: Local destination path
    :type destination: Path
    :param preserve_folder: If True, preserve the folder structure when downloading
    :type preserve_folder: bool
    :return: Updated URI and destination path
    :rtype: tuple
    """
    # List contents of S3 URI
    dir_contents = subprocess.run(["aws", "s3", "ls", uri], capture_output=True, text=True, check=True)
    
    entries = []

    # Iterate through the output lines and extract file/folder names
    for line in dir_contents.stdout.splitlines():
        parts = line.split()
        if len(parts) == 0:
            continue
    
        if parts[0] == "PRE":
            entries.append(parts[1])
    
        elif len(parts) >= 4:
            entries.append(parts[-1])
    
    if not entries:
        raise ValueError(f"  No files found in {uri}")

    # Find most recent file/folder and change URI
    latest = sorted(entries)[-1]
    
    if not uri.endswith("/"):
        uri += "/"
    
    uri += latest

    # If preserve_folder is True, append the last part of the URI to the destination path
    if preserve_folder:
        destination = destination / uri.rstrip("/").split("/")[-1]

    return uri, destination

def aws_download(uri, destination, include=None, recursive=False, most_recent=False, preserve_folder=False, dry_run=False, **kwargs):
    """
    Download files from AWS S3.

    :param uri: S3 URI to download from.
    :type uri: str
    :param destination: Local destination path.
    :type destination: Path
    :param include: Pattern to include specific files.
    :type include: str
    :param recursive: If True, download recursively.
    :type recursive: bool
    :param most_recent: If True, download only the most recent file.
    :type most_recent: bool
    :param preserve_folder: If True, preserve the folder structure when downloading.
    :type preserve_folder: bool
    :param dry_run: If True, show what would be downloaded without downloading files.
    :type dry_run: bool
    """
    destination = Path(destination)

    # Find most recent file/folder if requested
    if most_recent:
        uri, destination = find_latest_by_name(uri, destination, preserve_folder=preserve_folder)

    if dry_run:
        print(f"  Would download from {uri} to {destination}")
        return
    else:
        print(f"  Downloading from {uri} to {destination}")

        # Copy files from S3
        cmd = ["aws", "s3", "cp", uri, str(destination)]

        if recursive:
            cmd.append("--recursive")

        if include:
            cmd.extend([
                "--exclude", "*",
                "--include", include
            ])

        try:
            subprocess.run(cmd, check=True, text=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Download failed: {e.stderr}")
        except FileNotFoundError:
            raise RuntimeError("AWS CLI is not installed or not in PATH")
        except Exception as e:
            raise RuntimeError(f"Downloading failed: {e}")

def download_packages(packages, dry_run=False):
    """
    Download required packages for python and from a github repo

    :param packages: Packages section of YAML config
    :type packages: dictionary
    :param dry_run: If True, show what would be downloaded without downloading files.
    :type dry_run: bool
    """

    for package_type, package_list in packages.items():

        # If packages are for python pip install packages
        if package_type == "python":
            for package in package_list:
                if dry_run:
                    print(f"  Would have downloaded {package}")
                else:
                    try:
                        subprocess.run(["pip", "install", package], check=True)
                        print(f"  Downloaded {package}")
                    except subprocess.CalledProcessError:
                        raise RuntimeError("Ensure you have a python venv activated")

        # If packages are from git, find and install
        if package_type == "git":
            for package in package_list:
                download_git_package(**package, dry_run=dry_run)

def download_git_package(provider, owner, repo, destination=None, most_recent=True, version=None, include=None, unzip=False, dry_run=False):
    """
    Downloads git package that has tags. Will either grabs latest release or version given.
    Doesnt currently only works for public github repos

    :param provider: Either github or bitbucket depending on git host
    :type provider: str
    :param owner: Owner of repo eg GeoscienceAustralia
    :type owner: str
    :param repo: Name of repo eg dynadjust
    :type repo: str
    :param destination: Path when package should be downloaded
    :type destination: str or Path
    :param most_recent: If True, download only the most recent release.
    :type most_recent: bool
    :param version: Name of release to be downloaded
    :type version: str
    :param include: Pattern to include specific files.
    :type include: str
    :param unzip: Wether to unzip downloaded files
    :type unzip: bool
    :param dry_run: If True, show what would be downloaded without downloading files.
    :type dry_run: bool    
    """
    # Create URL depending on provider
    if provider == "github":

        url = f"https://api.github.com/repos/{owner}/{repo}/releases/"
        assets = "assets"
        if most_recent:
            url = url + "latest"
        elif version is not None:
            url = url + str(version)

    elif provider == "bitbucket":

        url = f"https://api.bitbucket.org/2.0/repositories/{owner}/{repo}/downloads"
        assets = "values"
        print( f"Bitbucket is currently not supported! Please manually download {repo}")
        return

    else:
        raise ValueError(f"Provider must be either github or bitbucket for {repo}")

    # Get API results
    release = requests.get(url)
    release.raise_for_status()

    files = []

    # Find all files to be downloaded
    for asset in release.json()[assets]:
        if include is None or include in asset["name"]:
            files.append(asset)

    if len(files) == 0:
        raise ValueError(
            f"No files found to download at {url}"
        )

    # Iterate through each file to be downloaded
    for file in files:
        output_path = Path(destination) / file["name"] if destination else Path(file["name"])

        # Check if file already exists
        if unzip:
            output_path_unzipped = output_path.with_suffix("")
            if output_path_unzipped.is_dir():
                print(f"  {output_path_unzipped} already exists")
                continue               
        else:
            if output_path.is_file():
                print(f"  {output_path} already exists")
                continue

        # Download files unless dry run is set
        if dry_run:
            print(f"  Would have downloaded package at {file["browser_download_url"]}")
        else:
            
            with requests.get(file["browser_download_url"], stream=True) as r:
                r.raise_for_status()

                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            print(f"  Downloaded package at {file["browser_download_url"]}")

            # If needed unzip file
            if unzip:
                try:
                    print(f"  Unzipping {output_path}")
                    shutil.unpack_archive(output_path)
                    print(f"  Unzipped {output_path}")
                    Path.unlink(output_path)
                except:
                    raise RuntimeError(f"Failed to unzip {output_path}")

    return output_path


def confirm_clean(processes):
    """
    Confirm with the user before cleaning specified processes.

    :param processes: List of processes to clean.
    :type processes: list
    """

    print("\nWARNING")
    print("The following processes will be permanently cleaned:")

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
            folder_path = (repo_root / folder).resolve()

            if not folder_path.is_relative_to(repo_root):
                raise ValueError(f"  Will not delete file outside repo root: {folder_path}")

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

if __name__ == "__main__":
    main()