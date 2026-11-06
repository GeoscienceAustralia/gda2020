from dataclasses import dataclass
from datetime import datetime, timezone
import argparse

parser = argparse.ArgumentParser(
        description="Dispatch batch jobs when AUSPOS downloads are finished."
    )

parser.add_argument(
    "--run-id",
    required=True,
    help="Unique identifier shared by the complete NGCA run.",
)

@dataclass(frozen=True)
class DispatchResult:
    jurisdiction: str
    submitted: bool
    job_id: str | None
    message: str

def dispatch_postprocessing(juris,archive,run_id,):
    """
    Dispatch post-processing for one completed jurisdiction.

    This is a placeholder. It will later call AWS Batch and return
    the submitted Batch job ID.
    """
    juris = juris.lower()

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