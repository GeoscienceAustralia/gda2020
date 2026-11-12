#!/usr/bin/env python3

"""Package NGCA jurisdiction data and optionally dispatch post-processing."""

import json
import os
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile
import sys

import boto3

from ngca_path import NGCA_DIR


@dataclass(frozen=True)
class DispatchResult:
    jurisdiction: str
    uploaded: bool
    submitted: bool
    job_id: str | None
    s3_key: str
    message: str


def package_filename(juris, run_id, partial=False):
    """Return the normal or partial package filename."""
    juris = juris.lower()
    run_suffix = run_id.removeprefix("ngca-")

    partial_suffix = "-partial" if partial else ""

    return (f"ngca-{juris}-{run_suffix}{partial_suffix}.zip")


def zip_and_upload_to_s3(folders, bucket, s3_key, zip_path):
    """ZIP several folders and upload the archive to S3."""
    zip_file = Path(zip_path)

    zip_file.parent.mkdir(parents=True, exist_ok=True)

    zip_file.unlink(missing_ok=True)

    uploaded = False

    # Try to create zip file containing all files in folders
    try:
        with ZipFile(zip_file, "w", compression=ZIP_DEFLATED, allowZip64=True) as archive:
            for folder_name in folders:
                folder = Path(folder_name)

                if not folder.is_dir():
                    raise FileNotFoundError(f"Folder does not exist: {folder}")

                for file in sorted(folder.rglob("*")):
                    if file.is_file():
                        archive.write(file, arcname=(Path(folder.name) / file.relative_to(folder)),)

        if not zip_file.is_file():
            raise RuntimeError(f"ZIP was not created: {zip_file}")

        if zip_file.stat().st_size == 0:
            raise RuntimeError(f"ZIP is empty: {zip_file}")

        print(f"Created package: {zip_file} ({zip_file.stat().st_size:,} bytes)", flush=True)

        # Upload zip file to s3
        boto3.client("s3").upload_file(str(zip_file), bucket, s3_key)

        uploaded = True

        print(f"Uploaded package: s3://{bucket}/{s3_key}", flush=True)

    finally:
        if uploaded:
            # Remove zip file if uploaded
            zip_file.unlink(missing_ok=True)

            print(f"Removed local package: {zip_file}", flush=True)

        elif zip_file.exists():
            print(f"Package retained after failure: {zip_file}", file=sys.stderr, flush=True)

def submit_postprocessing_job(juris, run_id):
    """Submit one post-processing job to AWS Batch."""
    # Get environment paramters for queue name and job definition
    queue = os.environ.get("POSTPROCESS_QUEUE", "").strip()
    job_definition = os.environ.get("POSTPROCESS_JOB_DEFINITION", "").strip()

    if not queue:
        raise RuntimeError("POSTPROCESS_QUEUE is not configured")

    if not job_definition:
        raise RuntimeError("POSTPROCESS_JOB_DEFINITION is not configured")

    job_name = f"{run_id}-{juris.lower()}-post"

    # Send postprocessing job
    response = boto3.client("batch").submit_job(
        jobName=job_name,
        jobQueue=queue,
        jobDefinition=job_definition,
        parameters={"jurisdiction": juris.lower(), "runId": run_id,},
        tags={
            "Workflow": "NGCA",
            "RunId": run_id,
            "Jurisdiction": juris.upper(),
            "Stage": "Postprocessing",
        },
    )

    job_id = str(response.get("jobId") or "").strip()

    if not job_id:
        raise RuntimeError(f"AWS Batch accepted the request but returned no job ID: {json.dumps(response)}")

    print(f"Submitted post-processing job: name={job_name} job_id={job_id}", flush=True)

    return job_id


def dispatch_postprocessing(juris, archive, run_id, s3_bucket, submit_job=False, partial=False):
    """Package, upload and optionally submit one jurisdiction.

    A partial package is uploaded for recovery only. Partial packages
    are never submitted to post-processing.
    """
    juris = juris.lower()

    if not s3_bucket:
        raise ValueError("s3_bucket is required")

    if not run_id:
        raise ValueError("run_id is required")

    if partial and submit_job:
        raise ValueError("A partial package cannot be submitted to post-processing")

    juris_dir = NGCA_DIR / juris
    sinex_dir = juris_dir / "sinexFiles"
    rinexantls_dir = juris_dir / "rinexantls"
    archive_dir = juris_dir / archive
    package_dir = juris_dir / "packages"

    # Get name of zip file
    zip_name = package_filename(juris, run_id, partial)

    zip_path = package_dir / zip_name

    s3_key = (f"NGCA/batch_jobs/{juris}/{run_id}/{zip_name}")

    folders_to_zip = [archive_dir, rinexantls_dir, sinex_dir]

    print()
    print("=" * 70)
    print(f"Packaging {'partial ' if partial else ''}post-processing input for {juris.upper()}")

    # Create zip file and upload to s3 if bucket given
    zip_and_upload_to_s3(folders_to_zip, s3_bucket, s3_key, zip_path)

    submitted = False
    job_id = None

    # Submit prost processing
    if submit_job:
        print(f"Dispatching post-processing for {juris.upper()}", flush=True)

        job_id = submit_postprocessing_job(juris, run_id)

        submitted = True

    print("=" * 70)

    if partial:
        message = ("Partial package uploaded; post-processing was not submitted")

    elif submitted:
        message = ("Package uploaded and post-processing job submitted")

    else:
        message = ("Package uploaded; post-processing dispatch was not requested")

    return DispatchResult(
        jurisdiction=juris,
        uploaded=True,
        submitted=submitted,
        job_id=job_id,
        s3_key=s3_key,
        message=message,
    )