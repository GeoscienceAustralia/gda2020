from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import boto3
from ngca_path import NGCA_DIR

BUCKET = "gda2020-long-term-archive"

@dataclass(frozen=True)
class DispatchResult:
    jurisdiction: str
    submitted: bool
    job_id: str | None
    message: str

def dispatch_postprocessing(juris, archive, run_id, s3_bucket):
    """
    Dispatch post-processing for one completed jurisdiction.

    This is a placeholder. It will later call AWS Batch and return
    the submitted Batch job ID.
    """
    run_id_sections = run_id.split("-")

    juris_run_id = f"{run_id_sections[0]}-{juris.upper()}-{run_id_sections[1]}"
    zip_name = f"{juris_run_id}.zip"
    
    juris = juris.lower()
    juris_dir = NGCA_DIR / juris
    sinex_dir = juris_dir / 'sinexFiles'
    rinexantls_dir = juris_dir / 'rinexantls'
    archive_dir = juris_dir / archive
    backup_dir = juris_dir / 'backup'
    zip_path = backup_dir / zip_name


    folders_to_zip = [sinex_dir, rinexantls_dir, archive_dir]
    s3_key = f"NGCA/batch_jobs/{juris}/${run_id}/{zip_name}"

    zip_and_upload_to_s3(folders_to_zip, s3_bucket, s3_key, zip_path)

    print()
    print("=" * 70)
    print(f"{juris.upper()} IS READY FOR POST-PROCESSING")
    print(f"Run ID:  {run_id}")
    print(f"Archive: {archive}")
    print(
        f"FUTURE ACTION: submit post-processing Batch job "
        f"for {juris.upper()}"
    )
    print("=" * 70)

    return DispatchResult(
        jurisdiction=juris,
        submitted=True,
        job_id=None,
        message="Placeholder dispatcher completed",
    )

def zip_and_upload_to_s3(folders, bucket, s3_key, zip_path):
    """
    Zip multiple folders and upload the resulting archive to S3.

    Example:
        zip_and_upload_to_s3(
            folders=["./data", "./config", "./logs"],
            bucket="my-backup-bucket",
            s3_key="backups/package.zip",
        )
    """

    zip_file = Path(zip_path)

    # Create the ZIP archive
    with ZipFile(zip_file, "w", compression=ZIP_DEFLATED) as archive:
        for folder_name in folders:
            folder = Path(folder_name)

            if not folder.is_dir():
                raise ValueError(f"Folder does not exist: {folder}")

            for file in folder.rglob("*"):
                if file.is_file():
                    # Produces paths such as:
                    # data/file.csv
                    # config/settings.json
                    archive.write(
                        file,
                        arcname=Path(folder.name) / file.relative_to(folder),
                    )

    # Upload to S3
    s3 = boto3.client("s3")

    try:
        s3.upload_file(
            str(zip_file),
            bucket,
            s3_key,
        )
    finally:
        # Don't leave the temporary ZIP behind
        zip_file.unlink(missing_ok=True)