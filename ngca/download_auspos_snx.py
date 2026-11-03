#!/usr/bin/env python3
"""Poll for and download AUSPOS SINEX files in parallel.

Reads successful job numbers from each jurisdiction's auspos_jobs.jsonl.
For job J, checks and downloads:
    https://d1hjwnsodwzynv.cloudfront.net/uploads/J/J.SNX

Files are written to:
    <root>/<juris>/<archive>/sinexfolder/J.SNX
"""

import argparse
import concurrent.futures
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests

BASE_URL = "https://d1hjwnsodwzynv.cloudfront.net/uploads"
DEFAULT_JURIS = ("tas", "act", "vic", "sa", "nt", "wa", "nsw", "qld")


@dataclass(frozen=True)
class Job:
    """
    Represents each job submitted to AUSPOS
    """
    juris: str
    number: str
    cluster_name: str
    destination: Path

    @property
    def url(self):
        """
        The URL that the sinex file can be retrived from
        """
        return f"{BASE_URL}/{self.number}/{self.number}.SNX"

    @property
    def output(self):
        """
        The destination output of the sinex file including renaming to the cluster ID name
        """
        name = self.cluster_name

        if name.endswith("_ls"):
            name = name[:-3]

        return self.destination / f"{name}.SNX"


def arguments():
    parser = argparse.ArgumentParser(
        description="Poll CloudFront and download completed AUSPOS SINEX files."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--archive", required=True)
    parser.add_argument("-j", "--juris", nargs="+", default=list(DEFAULT_JURIS))
    parser.add_argument("--initial-wait-minutes", type=float, default=30.0)
    parser.add_argument("--poll-minutes", type=float, default=10.0)
    parser.add_argument("--max-wait-hours", type=float, default=24.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--connect-timeout", type=float, default=20.0)
    parser.add_argument("--read-timeout", type=float, default=300.0)
    parser.add_argument(
        "--once", action="store_true",
        help="Check once immediately, download ready files, then exit."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Download again even if the final SNX file already exists."
    )
    return parser.parse_args()


def read_jobs(root, archive, jurisdictions, overwrite):
    """
    Create a list of jobs that need to have sinex downloaded. If overwrite include all jobs.
    Retruns dictionary of jobs.

    :param root: The root path of NGCA
    :type root: pathlib.Path
    :param archive:
    :type archive: str or pathlib.Path
    :param jurisdictions: Juris to be downloaded
    :type jurisdictions: List
    :param overwrite: Overwrite sinex files already downloaded
    :type overwrite bool:
    :return: jobs that need to be downloaded
    :rtype: Job
    """
    
    jobs = []
    missing_logs = 0

    # iterate through each juris
    for raw_juris in jurisdictions:
        juris = raw_juris.lower()
        archive_dir = root / juris / archive
        log = archive_dir / "auspos_jobs.jsonl"
        destination = root / juris / "sinexFiles"

        # Find .jsonl file with all job ids to be downloaded
        if not log.is_file():
            print(f"WARNING: no job log for {juris.upper()}: {log}", file=sys.stderr)
            missing_logs += 1
            continue

        # Build sinexFiles folder
        destination.mkdir(parents=True, exist_ok=True)

        seen = set()
        
        # Load each line of the json file, getting the number and cluster name of each job
        with log.open(encoding="utf-8") as source:
            for line_number, raw in enumerate(source, 1):
                try:
                    record = json.loads(raw)
                except json.JSONDecodeError as error:
                    print(f"WARNING: invalid JSON at {log}:{line_number}: {error}", file=sys.stderr)
                    continue

                number = str(record.get("job_number") or "").strip()
                cluster_name = str(record.get("cluster_name") or "").strip()

                if not record.get("ok") or not number or number in seen:
                    continue

                seen.add(number)

                job = Job(juris, number, cluster_name, destination)
                if overwrite or not job.output.is_file():
                    jobs.append(job)

    return jobs, missing_logs


def probe(job, connect_timeout: float):
    """
    Test if a job is ready to be downloaded or not. Looks for a status code of 200

    Returns (job, ready, status/error).

    :param job: job to be checked
    :type job: Job
    :param connect_timeout: Time until it gives up reaching url
    :type connect_timeout: float
    """
    try:
        response = requests.head(
            job.url,
            allow_redirects=True,
            timeout=(connect_timeout, connect_timeout),
            headers={"Cache-Control": "no-cache"},
        )
        if response.status_code == 200:
            return job, True, "HTTP 200"
        if response.status_code in (403, 404):
            return job, False, f"HTTP {response.status_code}"
        if response.status_code == 405:
            return job, None, "HEAD not supported (HTTP 405)"
        return job, False, f"HTTP {response.status_code}"
    except requests.RequestException as error:
        return job, False, f"probe error: {error}"


def download(job, connect_timeout, read_timeout):
    """
    Download the sinex file for a job

    Returns (job, ready, status/error).

    :param job: job to be checked
    :type job: Job
    :param connect_timeout: Time until it gives up reaching url
    :type connect_timeout: float
    :param read_timeout: Time after which to give up downloading
    :type read_timeout: float
    """
    # Create temp location to download sinex file to
    temp = job.output.with_suffix(job.output.suffix + ".part")

    # Try to download file
    try:
        temp.unlink(missing_ok=True)
        with requests.get(
            job.url,
            stream=True,
            allow_redirects=True,
            timeout=(connect_timeout, read_timeout),
            headers={"Cache-Control": "no-cache"},
        ) as response:
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "").lower()
            if "text/html" in content_type:
                raise RuntimeError(f"unexpected Content-Type: {content_type}")

            with temp.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        output.write(chunk)

        if temp.stat().st_size == 0:
            raise RuntimeError("downloaded file is empty")

        # Rename job
        temp.replace(job.output)
        return job, True, f"saved {job.output} ({job.output.stat().st_size:,} bytes)"
    
    except Exception as error:
        temp.unlink(missing_ok=True)
        return job, False, str(error)


def parallel_probe(jobs, workers, connect_timeout):
    """
    Check what jobs can be downloaded while paralleling

    :param jobs: List of jobs to be checked
    :type jobs: List of Job
    :param workers: number of parallel processes
    :type workers: int
    :param connect_timeout: Time after which to give up connecting
    :type connect_timeout: float
    """
    results = []
    
    # Create concurrent process and run probe in each
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(probe, job, connect_timeout) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    return results


def parallel_download(jobs, workers, connect_timeout, read_timeout):
    """
    Download sinex files for ready jobs in parallel

    :param jobs: List of jobs to be checked
    :type jobs: List of Job
    :param workers: number of parallel processes
    :type workers: int
    :param connect_timeout: Time after which to give up connecting
    :type connect_timeout: float
    :param read_timeout: Time after which to give up downloading
    :type read_timeout: float
    """
    results = []
    
    # Create concurrent process and run download in each
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(download, job, connect_timeout, read_timeout) for job in jobs
        ]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    return results

def check_completed_juris(pending, active_juris, completed_juris):
    """
    Check whether each jurisdiction has any jobs left in pending.

    If not, mark the jurisdiction complete and run juris_complete().

    :param pending: List of jobs that still need to be downloaded
    :type pending: List of Job
    :param juris: List of Juris being tested
    :type juris: List
    :param completed_juris: Set of jruis already completed
    :type completed_juris: Set
    """
    # Iterate through each juris
    for juris in active_juris:
        juris = juris.lower()

        if juris in completed_juris:
            continue

        # Find if any jobs in pending have a the juris
        if any(job.juris == juris for job in pending.values()):
            continue

        # If no jobs where found add juris to complete and run function
        completed_juris.add(juris)
        juris_complete(juris)

def juris_complete(juris):

    print(f"\n>>> {juris.lower()} is completely downloaded <<<")

def main():
    args = arguments()
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    if args.poll_minutes <= 0 and not args.once:
        raise SystemExit("--poll-minutes must be greater than zero")
    if args.max_wait_hours <= 0 and not args.once:
        raise SystemExit("--max-wait-hours must be greater than zero")

    # Find all jobs that need to be downloaded
    jobs, missing_logs = read_jobs(
        args.root.expanduser().resolve(), args.archive, args.juris, args.overwrite
    )
    
    # List of jobs still pending
    pending = {job.number: job for job in jobs}

    completed_juris = set()
    active_juris = {job.juris for job in jobs}

    print(f"Pending SINEX files: {len(pending)}; parallel workers: {args.workers}")
    if not pending:
        return 1 if missing_logs else 0

    # Start timer for total wait
    started = time.monotonic()
    deadline = started + args.max_wait_hours * 3600

    # Set intial wait time and sleep for that
    if not args.once and args.initial_wait_minutes > 0:
        seconds = args.initial_wait_minutes * 60
        print(f"Initial wait: {args.initial_wait_minutes:g} minute(s)")
        time.sleep(seconds)

    failed_downloads = set()


    while pending:
        now_text = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"\n[{now_text}] Checking {len(pending)} job(s)...")
        
        # Check all pending job to see which can be downloaded
        checks = parallel_probe(list(pending.values()), args.workers, args.connect_timeout)

        ready = []

        # Iterate through each job
        for job, is_ready, detail in sorted(checks, key=lambda item: (item[0].juris, item[0].number)):
            if is_ready:
                ready.append(job)
                #print(f"READY   {job.juris.upper()} {job.number}: {detail}")
            elif is_ready is None:
                print(f"ERROR   {job.juris.upper()} {job.number}: {detail}", file=sys.stderr)
            else:
                print(f"PENDING {job.juris.upper()} {job.number}: {detail}")

        # If there are files that are ready download them in parallel
        if ready:
            print(f"Downloading {len(ready)} ready file(s) in parallel...")
            results = parallel_download(
                ready, args.workers, args.connect_timeout, args.read_timeout
            )
            
            for job, ok, detail in sorted(results, key=lambda item: (item[0].juris, item[0].number)):
                if ok:
                    pending.pop(job.number, None)
                    #print(f"DONE    {job.juris.upper()} {job.number}: {detail}")
                else:
                    failed_downloads.add(job.number)
                    print(f"FAILED  {job.juris.upper()} {job.number}: {detail}", file=sys.stderr)

            check_completed_juris(pending, active_juris, completed_juris)

        # If there are no files left in pending then all have been downloaded
        if not pending:
            print("\nAll SINEX files have been downloaded.")
            return 0

        if pending:
            print(f"Still {len(pending)} job(s) left to be downloaded.")

        # If set to only run once then end
        if args.once:
            print(f"\nNot ready: {len(pending)} job(s).")
            return 2

        # If max wait time is reached then end
        if time.monotonic() >= deadline:
            print(f"\nMaximum wait reached; still pending: {len(pending)} job(s).", file=sys.stderr)
            for job in sorted(pending.values(), key=lambda item: (item.juris, item.number)):
                print(f"  {job.juris.upper()} {job.number}", file=sys.stderr)
            return 2

        # Sleep for set amount of time to retry
        sleep_seconds = min(args.poll_minutes * 60, max(0, deadline - time.monotonic()))
        print(f"Next check in {sleep_seconds / 60:g} minute(s).")
        time.sleep(sleep_seconds)

    return 1 if failed_downloads or missing_logs else 0


if __name__ == "__main__":
    raise SystemExit(main())
