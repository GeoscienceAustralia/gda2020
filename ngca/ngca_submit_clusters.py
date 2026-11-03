#!/usr/bin/env python3
"""Submit NGCA RINEX clusters in rinexantls/ to the AUSPOS API."""

import argparse
import json
import sys
import time
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

import requests

DEFAULT_API_URL = "https://gnss.ga.gov.au/ausposapi/jobs"
DEFAULT_JURIS = ("tas", "act", "vic", "sa", "nt", "wa", "nsw", "qld")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--archive", required=True)
    p.add_argument("-j", "--juris", nargs="+", default=list(DEFAULT_JURIS))
    p.add_argument("--email", required=True)
    p.add_argument("--api-url", default=DEFAULT_API_URL)
    p.add_argument("--timeout", type=int, default=1800)
    p.add_argument("--pause", type=float, default=0.1)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--stop-on-error", action="store_true")
    return p.parse_args()


def auspos_antenna_type(antenna, radome = "NONE"):
    """
    Create the 20-character ANT # / TYPE value expected by AUSPOS.

    RINEX antenna type occupies columns 1-16 and radome columns 17-20.
    Examples:
        LEIAX1202GG     NONE
        ASH701945C_M    NONE

    :param antenna: Name of antenna
    :type antenna: str
    :param radome: part of end of Antenna name "NONE"
    :type radome: str
    """
    antenna = antenna.strip()
    radome = (radome.strip() or "NONE")
    if len(antenna) > 16:
        raise ValueError(f"antenna name is longer than 16 characters: {antenna!r}")
    if len(radome) > 4:
        raise ValueError(f"radome code is longer than 4 characters: {radome!r}")
    return f"{antenna:<16}{radome:<4}"


def read_cluster(cluster_file, archive_dir):
    """
    Reads cluster file to find cluster name, antenna name, antenna height and rinex path

    :param cluster_file: cluster file e.g."16180_ls"
    :type cluster_file: pathlib.Path
    :param archive_dir: archive dir with rinex files
    :type archive_dir: pathlib.Path
    """
    entries = []
    with cluster_file.open(encoding="utf-8") as src:
        for line_number, raw in enumerate(src, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 3:
                raise ValueError(
                    f"{cluster_file}:{line_number}: expected filename, antenna type and height"
                )

            # ngca_prepFiles cluster format is: filename antenna height
            name, antenna, antenna_height = cols[0], cols[1], cols[2]
            
            # Format antenna_type str
            antenna_type = auspos_antenna_type(antenna, "NONE")

            # Find rinex path
            rinex_path = archive_dir / name
            if not rinex_path.is_file():
                raise FileNotFoundError(f"RINEX file not found: {rinex_path}")

            entries.append({
                "path": rinex_path,
                "name": name,
                "antenna_type": antenna_type,
                "antenna_height": antenna_height,
            })

    if not entries:
        raise ValueError(f"Cluster is empty: {cluster_file}")
    return entries


def successful_clusters(log_path):
    """
    Return set of clusters that have already been submitted to AUSPOS

    :param log_path: Path to log file
    :type log_path: pathlib.Path
    :return: set of clusters
    :rtype: set
    """
    completed = set()
    if not log_path.is_file():
        return completed
    with log_path.open(encoding="utf-8") as src:
        for raw in src:
            try:
                record = json.loads(raw)
            except json.JSONDecodeError:
                continue
            if record.get("ok") and record.get("cluster"):
                completed.add(record["cluster"])
    return completed


def append_record(jsonl_path: Path, record: dict):
    """
    :param jsonl_path: Path to jsonl file
    :type jsonl_path: pathlib.Path
    :param record: record to be appended
    :type record: dict
    """
    with jsonl_path.open("a", encoding="utf-8") as out:
        out.write(json.dumps(record, ensure_ascii=False) + "\n")

def submit_cluster(session, api_url, email, entries, timeout):
    """
    Submit a cluster to AUSPOS API and return http session and any response

    :param session: https session
    :type session: requests.Session()
    :param api_url: url to api
    :type api_url: str
    :param email: Email adress for AUSPOS responses
    :type email: str
    :param entries: list of files in cluster
    :type entries: List
    :param timeout: number of seconds to try http request for
    :type timeout: int
    """
    # Create the json playload for the cluster
    payload = {
        "files": [
            {"antenna_type": x["antenna_type"], "antenna_height": x["antenna_height"], "name": x["name"]}
            for x in entries
        ],
        "email": email,
    }

    # open required files
    with ExitStack() as stack:

        # Build multipart HTTP request
        multipart = {"payload": (None, json.dumps(payload), "application/json")}
        
        # Add each rinex file
        for item in entries:
            handle = stack.enter_context(item["path"].open("rb"))
            multipart[item["name"]] = (item["name"], handle, "text/plain")
        
        # Send HTTP request
        response = session.post(api_url, files=multipart, timeout=(30, timeout))

    try:
        body = response.json()
    except ValueError:
        body = {"raw_response": response.text}

    # Give good error message if failed
    if not response.ok:
        print("\n" + "=" * 70)
        print("AUSPOS SUBMISSION FAILED")
        print("=" * 70)
        print("HTTP status:", response.status_code)
        print("Response:")
        print(json.dumps(body, indent=2))
        print("\nPayload sent:")
        print(json.dumps(payload, indent=2))
        print("=" * 70)
        err = requests.HTTPError(
            f"{response.status_code} {response.reason}; AUSPOS response: {body}",
            response=response,
        )
        err.auspos_body = body
        raise err

    return response, body


def main():
    args = parse_args()
    root = args.root.expanduser().resolve()
    failures = 0

    # Create a HTTP session
    with requests.Session() as session:
        
        # Loop through each juris
        for juris_raw in args.juris:
            
            # Create required paths
            juris = juris_raw.lower()
            juris_dir = root / juris
            rinexantls_dir = juris_dir / "rinexantls"
            archive_dir = juris_dir / args.archive
            jsonl_path = archive_dir / "auspos_jobs.jsonl"

            # Check if rinexantls_dir exists
            if not rinexantls_dir.is_dir():
                print(f"ERROR: cluster directory not found: {rinexantls_dir}", file=sys.stderr)
                failures += 1
                continue

            # Find all cluster files
            clusters = sorted(
                p for p in rinexantls_dir.iterdir() 
                if p.is_file()
            )
            # Make list of files already submitted if using --resume
            done = successful_clusters(jsonl_path) if args.resume else set()
            
            print(f"{juris.upper()}: {len(clusters)} cluster file(s) found in {rinexantls_dir}")

            # Loop through each cluster
            for cluster_file in clusters:
                cluster_key = str(cluster_file.resolve())

                # skip cluster if already completed
                if cluster_key in done:
                    print(f"SKIP {juris.upper()} {cluster_file.name}: already submitted")
                    continue

                # Create job record
                record = {
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                    "juris": juris,
                    "archive": args.archive,
                    "cluster": cluster_key,
                    "cluster_name": cluster_file.name,
                    "ok": False,
                }

                try:
                    # Find all files in cluster
                    entries = read_cluster(cluster_file, archive_dir)

                    # Calculate size of cluster
                    total_bytes = sum(x["path"].stat().st_size for x in entries)
                    record["rinex_count"] = len(entries)
                    record["total_bytes"] = total_bytes
                    
                    #print(f"SUBMIT {juris.upper()} {cluster_file.name}: {len(entries)} RINEX, {total_bytes / 1024**3:.3f} GiB")
                    
                    # Dont submit to AUSPOS if dry run flag is set
                    if args.dry_run:
                        for x in entries:
                            print(f"  {x['name']}: antenna_type={x['antenna_type']!r}, height={x['antenna_height']}")
                        continue

                    # Submit to AUSPOS API
                    response, body = submit_cluster(session, args.api_url, args.email, entries, args.timeout)

                    # Get status and job number from response and save to json file
                    record["http_status"] = response.status_code
                    record["job_number"] = body.get("job_number") if isinstance(body, dict) else None
                    record["ok"] = True
                    append_record(jsonl_path, record)
                    
                    #print(f"  accepted: job_number={record.get('job_number')}")

                # If any error occurs also record it in json file
                except Exception as exc:
                    failures += 1
                    record["error"] = f"{type(exc).__name__}: {exc}"
                    if isinstance(exc, requests.HTTPError) and exc.response is not None:
                        record["http_status"] = exc.response.status_code
                        record["response"] = getattr(exc, "auspos_body", None)
                    if not args.dry_run:
                        append_record(jsonl_path, record)
                    print(f"  FAILED: {record['error']}", file=sys.stderr)
                    if args.stop_on_error:
                        return 1
                
                # Sleep for set time before sending next request
                if not args.dry_run and args.pause > 0:
                    time.sleep(args.pause)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
