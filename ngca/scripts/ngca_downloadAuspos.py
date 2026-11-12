#!/usr/bin/env python3

"""Poll for and download AUSPOS SINEX files.

Completed jurisdictions are packaged and uploaded when --s3-bucket is
provided. Post-processing is submitted only when
--dispatch-postprocessing is also provided.

At maximum wait, incomplete jurisdictions are uploaded as clearly named
partial packages and are never submitted to post-processing.
"""

import argparse
import concurrent.futures
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from ngca_path import NGCA_DIR


BASE_URL = (
    "https://d1hjwnsodwzynv.cloudfront.net/uploads"
)

DEFAULT_JURIS = (
    "tas",
    "act",
    "vic",
    "sa",
    "nt",
    "wa",
    "nsw",
    "qld",
)


@dataclass(frozen=True)
class Job:
    """One successful AUSPOS submission."""

    juris: str
    number: str
    cluster_name: str
    destination: Path

    @property
    def url(self):
        """The URL that the sinex file can be retrived from"""
        return (f"{BASE_URL}/{self.number}/{self.number}.SNX")

    @property
    def output(self):
        """The destination output of the sinex file including renaming to the cluster ID name"""
        name = self.cluster_name

        if name.endswith("_ls"):
            name = name[:-3]

        return self.destination / f"{name}.SNX"


def arguments():
    parser = argparse.ArgumentParser(description=("Poll CloudFront and download completed AUSPOS SINEX files."))

    parser.add_argument("--archive", required=True)
    parser.add_argument("-j", "--juris", nargs="+", default=list(DEFAULT_JURIS), choices=DEFAULT_JURIS, type=str.lower)
    parser.add_argument("--initial-wait-minutes", type=float, default=30.0)
    parser.add_argument("--poll-minutes", type=float, default=10.0)
    parser.add_argument("--max-wait-hours", type=float, default=24.0,)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--connect-timeout", type=float, default=20.0)
    parser.add_argument("--read-timeout", type=float, default=300.0)
    parser.add_argument("--run-id", help="Unique NGCA run identifier.")
    parser.add_argument("--s3-bucket", help=("S3 bucket for completed and partial post-processing packages."))
    parser.add_argument("--dispatch-postprocessing", action="store_true", help=("Submit the post-processing Batch job after a completed package is uploaded."))
    parser.add_argument("--once", action="store_true", help=("Check once immediately, download ready files, then exit."))
    parser.add_argument("--overwrite", action="store_true", help=("Download again even when the final SINEX file exists."))

    return parser.parse_args()


def validate_arguments(args):
    """Validate downloader configuration."""
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")

    if (args.poll_minutes <= 0 and not args.once):
        raise SystemExit("--poll-minutes must be greater than zero")

    if (args.max_wait_hours <= 0 and not args.once):
        raise SystemExit("--max-wait-hours must be greater than zero")

    if args.s3_bucket and not args.run_id:
        raise SystemExit("--run-id is required when --s3-bucket is used")

    if (args.dispatch_postprocessing and not args.s3_bucket):
        raise SystemExit("--s3-bucket is required when --dispatch-postprocessing is used")


def read_jobs(archive,jurisdictions,overwrite,):
    """ Create a list of jobs that need to have sinex downloaded. If overwrite include all jobs.
    Retruns dictionary of jobs."""
    jobs = []
    missing_logs = set()
    job_counts = {}

    # Iterate through each Juris
    for raw_juris in jurisdictions:
        juris = raw_juris.lower()

        juris_dir = NGCA_DIR / juris
        archive_dir = juris_dir / archive
        log_path = (archive_dir / "auspos_jobs.jsonl")
        destination = (juris_dir / "sinexFiles")

        # Find .jsonl file with all jobs to be downloaded
        if not log_path.is_file():
            print(f"WARNING: no job log for {juris.upper()}: {log_path}", file=sys.stderr)
            missing_logs.add(juris)
            continue

        # Build sinexFiles folder
        destination.mkdir(parents=True, exist_ok=True,)

        seen = set()
        job_counts[juris] = 0

        # Load each line of the jsonl file, getting the number and cluster name of each job
        with log_path.open(encoding="utf-8") as source:
            for line_number, raw in enumerate(source, 1):
                try:
                    record = json.loads(raw)

                except json.JSONDecodeError as error:
                    print(f"WARNING: invalid JSON at {log_path}:{line_number}: {error}", file=sys.stderr)
                    continue

                number = str(record.get("job_number") or "").strip()

                cluster_name = str(record.get("cluster_name") or "").strip()

                if (not record.get("ok") or not number or not cluster_name or number in seen):
                    continue

                seen.add(number)
                job_counts[juris] += 1

                job = Job(juris, number, cluster_name, destination)

                if (overwrite or not job.output.is_file()):
                    jobs.append(job)

    return jobs, missing_logs, job_counts


def probe(job, connect_timeout):
    """Check whether the AUSPOS SINEX file is ready for download"""
    try:
        response = requests.head(
            job.url,
            allow_redirects=True,
            timeout=(connect_timeout, connect_timeout),
            headers={"Cache-Control": "no-cache"},
        )

        if response.status_code == 200:
            return job, True, "HTTP 200"
        if response.status_code in (403,404,):
            return (job, False, f"HTTP {response.status_code}")
        if response.status_code == 405:
            return (job, True, "HEAD unsupported; attempting GET")
        return (job, False, f"HTTP {response.status_code}")

    except requests.RequestException as error:
        return (job, False, f"probe error: {error}")


def download(job, connect_timeout, read_timeout):
    """Download one SINEX file."""
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
            content_type = response.headers.get("Content-Type", "",).lower()

            if "text/html" in content_type:
                raise RuntimeError(f"unexpected Content-Type: {content_type}")

            with temp.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        output.write(chunk)

        if temp.stat().st_size == 0:
            raise RuntimeError("downloaded file is empty")

        # Rename job from temp location
        temp.replace(job.output)
        return (job, True, f"saved {job.output}")

    except Exception as error:
        temp.unlink(missing_ok=True)
        return (job, False, str(error))


def parallel_probe(jobs, workers, connect_timeout):
    """Probe SINEX availability concurrently."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(lambda job: probe(job, connect_timeout), jobs))

def parallel_download(jobs, workers, connect_timeout, read_timeout):
    """Download ready SINEX files concurrently."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(lambda job: download(job, connect_timeout, read_timeout), jobs))

def finalise_jurisdiction(juris, archive, run_id, s3_bucket, dispatch_enabled):
    """Upload completed input and optionally dispatch it."""
    if not s3_bucket:
        return {
            "jurisdiction": juris,
            "uploaded": False,
            "submitted": False,
            "message": "Files remain local",
        }

    # Import only if needed
    from ngca_dispatchPost import dispatch_postprocessing

    result = dispatch_postprocessing(
        juris=juris,
        archive=archive,
        run_id=run_id,
        s3_bucket=s3_bucket,
        submit_job=dispatch_enabled,
        partial=False,
    )

    return {
        "jurisdiction": juris,
        "uploaded": result.uploaded,
        "submitted": result.submitted,
        "message": result.message,
    }


def check_completed_juris(pending, active_juris, completed_juris, completion_results, archive, run_id, s3_bucket, dispatch_enabled):
    """Handle jurisdictions that have no pending downloads."""

    # Check each juris that is active
    for juris in sorted(active_juris):

        # If juris is already complete skip
        if juris in completed_juris:
            continue

        # If juris still has jobs to be downloaded skip
        if any(job.juris == juris for job in pending.values()):
            continue

        # Try to finalise juris if all jobs are completed
        try:
            result = finalise_jurisdiction(
                juris=juris,
                archive=archive,
                run_id=run_id,
                s3_bucket=s3_bucket,
                dispatch_enabled=dispatch_enabled,
            )

        except Exception as error:
            result = {
                "jurisdiction": juris,
                "uploaded": False,
                "submitted": False,
                "message": (f"{type(error).__name__}: {error}"),
            }

            print(f"ERROR: finalisation failed for {juris.upper()}: {result['message']}", file=sys.stderr)

        completion_results[juris] = result
        completed_juris.add(juris)

        if result["submitted"]:
            print(f"DISPATCHED {juris.upper()}")

        elif result["uploaded"]:
            print(f"UPLOADED {juris.upper()}: dispatch not requested")

        elif result["message"] == "Files remain local":
            print(f"COMPLETE {juris.upper()}: files remain local")


def upload_partial_jurisdictions(pending, archive, run_id, s3_bucket):
    """Upload incomplete timeout packages without dispatching."""
    if not s3_bucket:
        return

    from ngca_dispatchPost import dispatch_postprocessing

    pending_juris = sorted({job.juris for job in pending.values()})

    for juris in pending_juris:
        try:
            result = dispatch_postprocessing(
                juris=juris,
                archive=archive,
                run_id=run_id,
                s3_bucket=s3_bucket,
                submit_job=False,
                partial=True,
            )

            print(f"PARTIAL UPLOAD {juris.upper()}: {result.message}", file=sys.stderr,)

        except Exception as error:
            print(f"ERROR: partial upload failed for {juris.upper()}: {type(error).__name__}: {error}", file=sys.stderr)


def completed_exit_code(completion_results, missing_logs, dispatch_enabled):
    """Return the final code after all downloads complete."""
    if missing_logs:
        return 1

    failed_uploads = [
        result
        for result in completion_results.values()
        if (result["message"] != "Files remain local" and not result["uploaded"])
    ]

    if failed_uploads:
        return 3

    if dispatch_enabled:
        failed_dispatches = [
            result
            for result in completion_results.values()
            if not result["submitted"]
        ]

        if failed_dispatches:
            return 3

    return 0


def main():
    args = arguments()
    validate_arguments(args)

    # Find all jobs to be downloaded
    jobs, missing_logs, job_counts = read_jobs(
        archive=args.archive,
        jurisdictions=args.juris,
        overwrite=args.overwrite,
    )

    # Find if there are no jobs for a juris
    empty_logs = {
        juris.lower()
        for juris in args.juris
        if (juris.lower() not in missing_logs and job_counts.get(juris.lower(),0) == 0)
    }

    if empty_logs:
        for juris in sorted(empty_logs):
            print(f"ERROR: no successful AUSPOS jobs for {juris.upper()}", file=sys.stderr)

        missing_logs.update(empty_logs)

    # Find which juris have active jobs
    active_juris = {
        juris.lower()
        for juris in args.juris
        if (juris.lower() not in missing_logs and job_counts.get(juris.lower(), 0) > 0)
    }

    # List of all jobs still to be downloaded
    pending = {job.number: job for job in jobs}

    completed_juris = set()
    completion_results = {}

    print(f"Pending SINEX files: {len(pending)}; parallel workers: {args.workers}")

    check_completed_juris(
        pending=pending,
        active_juris=active_juris,
        completed_juris=completed_juris,
        completion_results=completion_results,
        archive=args.archive,
        run_id=args.run_id,
        s3_bucket=args.s3_bucket,
        dispatch_enabled=args.dispatch_postprocessing,
    )

    # If no pending jobs left, complete process
    if not pending:
        print("\nAll available SINEX files have been downloaded.")

        return completed_exit_code(
            completion_results,
            missing_logs,
            args.dispatch_postprocessing,
        )

    deadline = (time.monotonic() + args.max_wait_hours * 3600)

    # Wait inital wait time
    if (not args.once and args.initial_wait_minutes > 0):
        print(f"Initial wait: {args.initial_wait_minutes:g} minute(s)")

        time.sleep(args.initial_wait_minutes * 60)

    # While there are still jobs to be downloaded
    while pending:
        now_text = time.strftime("%Y-%m-%d %H:%M:%S")

        print(f"\n[{now_text}] Checking {len(pending)} job(s)...")

        # Check all pending jobs to see which can be downloaded
        checks = parallel_probe(list(pending.values()), args.workers, args.connect_timeout)

        ready = [job for job, is_ready, _ in checks if is_ready]

        for job, is_ready, detail in checks:
            if not is_ready:
                print(f"PENDING {job.juris.upper()} {job.number}: {detail}")

        # If there are files that are available, download them in parallel
        if ready:
            print(f"Downloading {len(ready)} ready file(s)...")

            results = parallel_download(ready, args.workers, args.connect_timeout, args.read_timeout)

            for job, succeeded, detail in results:
                if succeeded:
                    # Remove job from pending if downloaded
                    pending.pop(job.number, None)

                else:
                    print(f"FAILED {job.juris.upper()} {job.number}: {detail}", file=sys.stderr)

            check_completed_juris(
                pending=pending,
                active_juris=active_juris,
                completed_juris=completed_juris,
                completion_results=(completion_results),
                archive=args.archive,
                run_id=args.run_id,
                s3_bucket=args.s3_bucket,
                dispatch_enabled=args.dispatch_postprocessing
            )

        # If there are no jobs left to be downloaded
        if not pending:
            print("\nAll SINEX files have been downloaded.")

            return completed_exit_code(completion_results, missing_logs, args.dispatch_postprocessing)

        if args.once:
            upload_partial_jurisdictions(pending, args.archive, args.run_id, args.s3_bucket,)
            return 2

        # If total wait time has been hit
        if time.monotonic() >= deadline:
            print(f"\nMaximum wait reached; still pending: {len(pending)} job(s).", file=sys.stderr)

            for job in sorted(pending.values(), key=lambda item: (item.juris, item.number),):
                print(f"  {job.juris.upper()} {job.number}", file=sys.stderr,)

            upload_partial_jurisdictions(pending, args.archive, args.run_id, args.s3_bucket)

            return 2

        sleep_seconds = min(args.poll_minutes * 60, max(0, deadline - time.monotonic()),)

        print(f"Still {len(pending)} job(s) left to download.")

        print(f"Next check in {sleep_seconds / 60:g} minute(s).")

        time.sleep(sleep_seconds)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())