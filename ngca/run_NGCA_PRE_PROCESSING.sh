#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob
export PYTHONUNBUFFERED=1

readonly SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1
    pwd
)"

readonly -a ALL_JURIS=(
    tas
    act
    vic
    sa
    nt
    wa
    nsw
    qld
)

JURIS_LIST=()
EMAIL=""
RUN_ID=""
ARCHIVE=""
PROCESS_ALL=false
DISPATCH_POSTPROCESSING=false


# These can be supplied through the AWS Batch job definition.
AUSPOS_SUBMISSION_DELAY_SECONDS="${AUSPOS_SUBMISSION_DELAY_SECONDS:-0.1}"
NGCA_RUN_BUCKET="${NGCA_RUN_BUCKET:-}"
POSTPROCESS_QUEUE="${POSTPROCESS_QUEUE:-}"
POSTPROCESS_JOB_DEFINITION="${POSTPROCESS_JOB_DEFINITION:-}"

# Per-jurisdiction processing state.
declare -A JURIS_STATUS
declare -A JURIS_FAILED_STAGE
declare -A JURIS_MESSAGE
declare -A JURIS_NUM_RINEX
declare -A JURIS_NUM_JOBS

SUBMITTED_JURIS=()

NOTES_FILE=""
SUMMARY_JSON=""
SUMMARY_TSV=""

usage() {
    cat <<EOF
Usage:
  $(basename "$0") --juris JURIS [JURIS ...] --email ADDRESS [options]
  $(basename "$0") --all --email ADDRESS [options]

Options:
  -j, --juris JURIS...   Jurisdictions to process
      --all              Process all jurisdictions
  -e, --email ADDRESS    Email address supplied to AUSPOS
  -r, --run-id VALUE     Unique run identifier
      --dispatch-postprocessing
                         Submit a post-processing AWS Batch job after
                         each jurisdiction package is uploaded to S3.
                         Without this flag, packages are uploaded only.
  -h, --help             Show this help

Environment variables:
  AUSPOS_SUBMISSION_DELAY_SECONDS
      Delay between AUSPOS submissions. Default: 0.1 seconds.

  NGCA_RUN_BUCKET
      Optional bucket used for resumable state and packaged inputs

  POSTPROCESS_QUEUE
      AWS Batch post-processing job queue. Required when 
      --dispatch-postprocessing is used

  POSTPROCESS_JOB_DEFINITION
      AWS Batch post-processing job definition. Required when
      --dispatch-postprocessing is used
EOF
}

# Function for consistent log files generation
log() {
    local level="$1"
    local jurisdiction="$2"
    local stage="$3"

    shift 3

    printf '%s level=%s run_id=%s jurisdiction=%s stage=%s message=%s\n' \
        "$(date --iso-8601=seconds)" \
        "$level" \
        "$RUN_ID" \
        "${jurisdiction^^}" \
        "$stage" \
        "$*"
}

# Function for consistent fail line generation
fail_controller() {
    local stage="$1"
    local message="$2"
    local exit_code="${3:-1}"

    printf '%s level=ERROR run_id=%s stage=%s message=%s\n' \
        "$(date --iso-8601=seconds)" \
        "$RUN_ID" \
        "$stage" \
        "$message" >&2

    exit "$exit_code"
}

# Check that juris entered is a valid juris
is_valid_jurisdiction() {
    local requested="$1"
    local valid

    for valid in "${ALL_JURIS[@]}"; do
        if [[ "$requested" == "$valid" ]]; then
            return 0
        fi
    done

    return 1
}

# Check if juris is valid, ignore duplicates and make all lowercase
add_jurisdiction() {
    local requested="${1,,}"
    local existing

    if ! is_valid_jurisdiction "$requested"; then
        printf 'ERROR: Unknown jurisdiction: %s\n' \
            "$requested" >&2
        exit 64
    fi

    # Ignore duplicated jurisdiction arguments.
    for existing in "${JURIS_LIST[@]}"; do
        if [[ "$existing" == "$requested" ]]; then
            return 0
        fi
    done

    JURIS_LIST+=("$requested")
}

# Change variables to show jurisdiction is processing
mark_jurisdiction_running() {
    local jurisdiction="$1"

    JURIS_STATUS["$jurisdiction"]="RUNNING"
    JURIS_FAILED_STAGE["$jurisdiction"]=""
    JURIS_MESSAGE["$jurisdiction"]=""

    log INFO "$jurisdiction" START \
        "Jurisdiction processing started"
}

# Change variables to show jurisdiction has failed
mark_jurisdiction_failed() {
    local jurisdiction="$1"
    local stage="$2"
    local message="$3"

    JURIS_STATUS["$jurisdiction"]="FAILED"
    JURIS_FAILED_STAGE["$jurisdiction"]="$stage"
    JURIS_MESSAGE["$jurisdiction"]="$message"

    log ERROR "$jurisdiction" "$stage" "$message"
}

# Change variables to show jurisdiction has successed
mark_jurisdiction_submitted() {
    local jurisdiction="$1"

    JURIS_STATUS["$jurisdiction"]="AUSPOS_SUBMITTED"
    JURIS_FAILED_STAGE["$jurisdiction"]=""
    JURIS_MESSAGE["$jurisdiction"]=""

    SUBMITTED_JURIS+=("$jurisdiction")
}

# Change variables to show jurisdiction has successed
mark_jurisdiction_complete() {
    local jurisdiction="$1"
    local message

    JURIS_STATUS["$jurisdiction"]="SUCCEEDED"
    JURIS_FAILED_STAGE["$jurisdiction"]=""
    JURIS_MESSAGE["$jurisdiction"]=""

    if "$DISPATCH_POSTPROCESSING"; then
        message="SINEX files downloaded, package uploaded, and post-processing job submitted"
    elif [[ -n "$NGCA_RUN_BUCKET" ]]; then
        message="SINEX files downloaded and post-processing package uploaded to S3"
    else
        message="SINEX files downloaded; outputs remain local"
    fi

    log INFO "$jurisdiction" COMPLETE "$message"
}

# Counts files in a directory
count_regular_files() {
    local directory="$1"

    find "$directory" \
        -maxdepth 1 \
        -type f \
        -printf '.' |
        wc -c
}

# Counts rinex files in directory
count_rinex_files() {
    local directory="$1"

    find "$directory" \
        -maxdepth 1 \
        -type f \
        \( -iname '*.*o' -o -iname '*.rnx' \) \
        -printf '.' |
        wc -c
}

auspos_state_key() {
    local jurisdiction="$1"

    printf '%s/%s/%s/auspos_jobs.jsonl' \
        "$RUN_ID" \
        "$jurisdiction" \
        "$ARCHIVE"
}

restore_auspos_state() {
    local jurisdiction="$1"
    local jsonl_path
    local object_key

    if [[ -z "$NGCA_RUN_BUCKET" ]]; then
        log INFO "$jurisdiction" AUSPOS_STATE \
            "No NGCA_RUN_BUCKET configured; using local submission state"

        return 0
    fi

    jsonl_path="$SCRIPT_DIR/$jurisdiction/$ARCHIVE/auspos_jobs.jsonl"
    object_key="NGCA/batch_jobs/${jurisdiction}/${RUN_ID}/auspos_jobs.jsonl"

    if aws s3api head-object \
        --bucket "$NGCA_RUN_BUCKET" \
        --key "$object_key" \
        >/dev/null 2>&1
    then
        if ! aws s3 cp \
            "s3://${NGCA_RUN_BUCKET}/${object_key}" \
            "$jsonl_path"
        then
            log ERROR "$jurisdiction" AUSPOS_STATE \
                "Existing submission state could not be downloaded"

            return 1
        fi

        log INFO "$jurisdiction" AUSPOS_STATE \
            "Restored submission state from S3"
    else
        log INFO "$jurisdiction" AUSPOS_STATE \
            "No S3 submission state found; using local state if present"
    fi
}

save_auspos_state() {
    local jurisdiction="$1"
    local jsonl_path
    local object_key

    if [[ -z "$NGCA_RUN_BUCKET" ]]; then
        log INFO "$jurisdiction" AUSPOS_STATE \
            "No NGCA_RUN_BUCKET configured; submission state remains local"

        return 0
    fi

    jsonl_path="$SCRIPT_DIR/$jurisdiction/$ARCHIVE/auspos_jobs.jsonl"
    object_key="NGCA/batch_jobs/${jurisdiction}/${RUN_ID}/auspos_jobs.jsonl"

    if [[ ! -s "$jsonl_path" ]]; then
        log ERROR "$jurisdiction" AUSPOS_STATE \
            "No non-empty submission log exists to upload"

        return 1
    fi

    if ! aws s3 cp \
        "$jsonl_path" \
        "s3://${NGCA_RUN_BUCKET}/${object_key}"
    then
        log ERROR "$jurisdiction" AUSPOS_STATE \
            "Submission state upload failed"

        return 1
    fi

    log INFO "$jurisdiction" AUSPOS_STATE \
        "Saved submission state to S3"
}

# Main function to process jurisdiction
process_jurisdiction() {
    local jurisdiction="$1"

    local archive_directory
    local antenna_directory
    local num_jobs
    local num_rinex
    local stage_exit

    archive_directory="$SCRIPT_DIR/$jurisdiction/$ARCHIVE"
    antenna_directory="$SCRIPT_DIR/$jurisdiction/rinexantls"

    mark_jurisdiction_running "$jurisdiction"

    #
    # Stage 1: Validate downloaded input directories
    #

    if [[ ! -d "$archive_directory" ]]; then
        mark_jurisdiction_failed "$jurisdiction" DOWNLOAD_VALIDATION "Archive directory not found: $archive_directory"

        return 1
    fi

    if [[ ! -d "$antenna_directory" ]]; then
        mark_jurisdiction_failed "$jurisdiction" DOWNLOAD_VALIDATION "Cluster directory not found: $antenna_directory"

        return 1
    fi

    #
    # Stage 2: Run verifySub
    #

    log INFO "$jurisdiction" VERIFY_SUB "Running ngca_verifySub.py"

    if (
        cd "$archive_directory"
        ngca_verifySub.py
    ); then
        log INFO "$jurisdiction" VERIFY_SUB "ngca_verifySub.py completed successfully"
    else
        stage_exit=$?

        mark_jurisdiction_failed "$jurisdiction" VERIFY_SUB "ngca_verifySub.py failed with exit code $stage_exit"

        return 1
    fi

    #
    # Stage 3: Run prepFiles
    #

    log INFO "$jurisdiction" PREP_FILES "Running ngca_prepFiles.py"

    if (
        cd "$archive_directory"
        ngca_prepFiles.py --juris "$jurisdiction"
    ); then
        log INFO "$jurisdiction" PREP_FILES "ngca_prepFiles.py completed successfully"
    else
        stage_exit=$?

        mark_jurisdiction_failed "$jurisdiction" PREP_FILES "ngca_prepFiles.py failed with exit code $stage_exit"

        return 1
    fi

    #
    # Stage 4: Validate preprocessing output
    #

    num_jobs="$(count_regular_files "$antenna_directory")"
    num_rinex="$(count_rinex_files "$archive_directory")"

    JURIS_NUM_JOBS["$jurisdiction"]="$num_jobs"
    JURIS_NUM_RINEX["$jurisdiction"]="$num_rinex"

    if ((num_rinex == 0)); then
        mark_jurisdiction_failed "$jurisdiction" PREP_VALIDATION "No RINEX files found in $archive_directory"

        return 1
    fi

    if ((num_jobs == 0)); then
        mark_jurisdiction_failed "$jurisdiction" PREP_VALIDATION "No AUSPOS cluster files found in $antenna_directory"

        return 1
    fi

    printf '%s has %d RINEX and %d jobs\n' "${jurisdiction^^}" "$num_rinex" "$num_jobs" >> "$NOTES_FILE"

    log INFO "$jurisdiction" PREP_COMPLETE "Preprocessing complete rinex=$num_rinex jobs=$num_jobs"

    #
    # Stage 5: Submit clusters to AUSPOS
    #

    log INFO "$jurisdiction" AUSPOS_SUBMISSION "Starting serial AUSPOS submission"

    if ! restore_auspos_state "$jurisdiction"; then
        mark_jurisdiction_failed "$jurisdiction" AUSPOS_STATE_RESTORE "Could not restore AUSPOS submission state"
        return 1
    fi

    submit_exit=0

    if ngca_submitClusters.py \
        --juris "$jurisdiction" \
        --archive "$ARCHIVE" \
        --email "$EMAIL" \
        --pause "$AUSPOS_SUBMISSION_DELAY_SECONDS" \
        --resume \
        #--dry-run # tmp while testing
    then
        submit_exit=0
    else
        submit_exit=$?

    fi

    # Save the complete or partial auspos submission list
    if ! save_auspos_state "$jurisdiction"; then
        mark_jurisdiction_failed "$jurisdiction" AUSPOS_STATE_UPLOAD "Could not save AUSPOS submission state"

        return 1
    fi

    if ((submit_exit != 0)); then
        mark_jurisdiction_failed "$jurisdiction" AUSPOS_SUBMISSION "submitClusters failed with exit code $submit_exit"

        return 1
    fi

    log INFO "$jurisdiction" AUSPOS_SUBMISSION "Every expected cluster has been submitted"

    mark_jurisdiction_submitted "$jurisdiction"
    return 0
}

run_auspos_downloader() {
    local downloader_exit
    local jurisdiction
    local -a downloader_args

    if ((${#SUBMITTED_JURIS[@]} == 0)); then
        log ERROR NONE AUSPOS_DOWNLOAD \
            "No jurisdictions were successfully submitted"

        return 1
    fi

    log INFO NONE AUSPOS_DOWNLOAD \
        "Monitoring jurisdictions: ${SUBMITTED_JURIS[*]^^}"

    downloader_args=(
        --archive "$ARCHIVE" \
        --juris "${SUBMITTED_JURIS[@]}" \
        --initial-wait-minutes \
            "${AUSPOS_INITIAL_WAIT_MINUTES:-1}" \
        --poll-minutes \
            "${AUSPOS_POLL_MINUTES:-10}" \
        --max-wait-hours \
            "${AUSPOS_MAX_WAIT_HOURS:-24}" \
        --workers \
            "${AUSPOS_DOWNLOAD_WORKERS:-2}" \
        --run-id "$RUN_ID" \
    )

    if [[ -n "$NGCA_RUN_BUCKET" ]];; then
        downloader_args+=( --s3-bucket "$NGCA_RUN_BUCKET")
    fi

    if "$DISPATCH_POSTPROCESSING"; then
        downloader_args+=( --dispatch-postprocessing)

        log INFO NONE POSTPROCESS_DISPATCH "Automatic post-processing enabled queue=$POSTPROCESS_QUEUE"
    else
        log INFO NONE POSTPROCESS_DISPATCH "Automatic post_processing dispatch disabled"
    fi

    if ngca_downloadAuspos.py "${downloader_args[@]}"; then
        downloader_exit=0

        for jurisdiction in "${SUBMITTED_JURIS[@]}"; do
            mark_jurisdiction_complete "$jurisdiction"
        done
    else
        downloader_exit=$?

        case "$downloader_exit" in
            1)
                log ERROR NONE AUSPOS_DOWNLOAD \
                    "One or more AUSPOS job logs were missing"
                ;;

            2)
                log ERROR NONE AUSPOS_DOWNLOAD \
                    "One or more AUSPOS jobs did not finish before the deadline"
                ;;

            3)
                log ERROR NONE POSTPROCESS_DISPATCH \
                    "One or more jurisdictions could not be dispatched to post-processing"
                ;;

            *)
                log ERROR NONE AUSPOS_DOWNLOAD \
                    "Downloader failed with unexpected exit code $downloader_exit"
                ;;
        esac

        # At this level, the Python downloader logs the specific
        # jurisdiction that failed. The shell only knows that the
        # combined downloader operation failed.
        for jurisdiction in "${SUBMITTED_JURIS[@]}"; do
            if [[ "${JURIS_STATUS[$jurisdiction]}" == "AUSPOS_SUBMITTED" ]]; then
                mark_jurisdiction_failed \
                    "$jurisdiction" \
                    AUSPOS_DOWNLOAD_OR_DISPATCH \
                    "Combined downloader failed with code $downloader_exit; check logs."
            fi
        done
    fi

    return "$downloader_exit"
}

# Print a summary to show which juris has passed or failed
print_summary() {
    local jurisdiction
    local status
    local passed=0
    local failed=0

    printf '\n'
    printf '%s\n' \
        '================================================================'
    printf 'NGCA CONTROLLER SUMMARY\n'
    printf 'Run ID:  %s\n' "$RUN_ID"
    printf 'Archive: %s\n' "$ARCHIVE"
    printf '%s\n' \
        '================================================================'

    # Iterate through each juris to check if successed or failed and reason
    for jurisdiction in "${JURIS_LIST[@]}"; do
        status="${JURIS_STATUS[$jurisdiction]:-UNKNOWN}"

        case "$status" in
            SUCCEEDED)
                ((passed += 1))

                printf '%-4s  %-10s  RINEX=%-6s  JOBS=%-6s\n' \
                    "${jurisdiction^^}" \
                    "SUCCEEDED" \
                    "${JURIS_NUM_RINEX[$jurisdiction]:-unknown}" \
                    "${JURIS_NUM_JOBS[$jurisdiction]:-unknown}"
                ;;

            FAILED)
                ((failed += 1))

                printf '%-4s  %-10s  STAGE=%-20s  %s\n' \
                    "${jurisdiction^^}" \
                    "FAILED" \
                    "${JURIS_FAILED_STAGE[$jurisdiction]:-UNKNOWN}" \
                    "${JURIS_MESSAGE[$jurisdiction]:-No message}"
                ;;

            *)
                ((failed += 1))

                printf '%-4s  %-10s  STAGE=%-20s  %s\n' \
                    "${jurisdiction^^}" \
                    "UNKNOWN" \
                    "UNKNOWN" \
                    "Jurisdiction did not reach a terminal status"
                ;;
        esac
    done

    printf '%s\n' \
        '================================================================'
    printf 'Passed: %d\n' "$passed"
    printf 'Failed: %d\n' "$failed"
    printf '%s\n' \
        '================================================================'
}

#
# Parse command-line arguments
#

while (($# > 0)); do
    case "$1" in
        -j|--juris)
            shift

            if (($# == 0)) || [[ "$1" == -* ]]; then
                printf 'ERROR: --juris requires at least one value\n' \
                    >&2
                exit 64
            fi

            while (($# > 0)) && [[ "$1" != -* ]]; do
                add_jurisdiction "$1"
                shift
            done
            ;;

        --all)
            PROCESS_ALL=true
            shift
            ;;

        -e|--email)
            if (($# < 2)); then
                printf 'ERROR: %s requires a value\n' \
                    "$1" >&2
                exit 64
            fi

            EMAIL="$2"
            shift 2
            ;;

        -r|--run-id)
            if (($# < 2)); then
                printf 'ERROR: %s requires a value\n' \
                    "$1" >&2
                exit 64
            fi

            RUN_ID="$2"
            shift 2
            ;;

        --dispatch-postprocessing)
            DISPATCH_POSTPROCESSING=true
            shift
            ;;
        
        -h|--help)
            usage
            exit 0
            ;;

        --)
            shift
            break
            ;;

        *)
            printf 'ERROR: Unknown option: %s\n' \
                "$1" >&2
            usage >&2
            exit 64
            ;;
    esac
done

#
# Validate controller arguments
#

if "$PROCESS_ALL" && ((${#JURIS_LIST[@]} > 0)); then
    printf 'ERROR: Use either --all or --juris, not both\n' \
        >&2
    exit 64
fi

if "$PROCESS_ALL"; then
    JURIS_LIST=("${ALL_JURIS[@]}")
elif ((${#JURIS_LIST[@]} == 0)); then
    printf 'ERROR: Specify --juris or --all\n' >&2
    exit 64
fi

if [[ -z "$EMAIL" ]]; then
    printf 'ERROR: --email is required\n' >&2
    exit 64
fi

if [[ -z "$RUN_ID" ]]; then
    RUN_ID="ngca-$(date -u +%Y%m%dT%H%M)"
fi

if "$DISPATCH_POSTPROCESSING"; then
    if [[ -z "$NGCA_RUN_BUCKET" ]]; then
        fail_controller \
            STARTUP_VALIDATION \
            "--dispatch-postprocessing requires NGCA_RUN_BUCKET" \
            64
    fi

    if [[ -z "$POSTPROCESS_QUEUE" ]]; then
        fail_controller \
            STARTUP_VALIDATION \
            "--dispatch-postprocessing requires POSTPROCESS_QUEUE" \
            64
    fi

    if [[ -z "$POSTPROCESS_JOB_DEFINITION" ]]; then
        fail_controller \
            STARTUP_VALIDATION \
            "--dispatch-postprocessing requires POSTPROCESS_JOB_DEFINITION" \
            64
    fi
fi

if [[ ! "$AUSPOS_SUBMISSION_DELAY_SECONDS" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    printf 'ERROR: Invalid AUSPOS submission delay: %s\n' \
        "$AUSPOS_SUBMISSION_DELAY_SECONDS" >&2
    exit 64
fi

#
# Configure command path
#

if [[ -d "$SCRIPT_DIR/scripts" ]] &&
   [[ ":$PATH:" != *":$SCRIPT_DIR/scripts:"* ]]; then
    export PATH="$PATH:$SCRIPT_DIR/scripts"
fi

#
# Validate dependencies before changing anything
#

for required_command in \
    ngca_get.py \
    ngca_verifySub.py \
    ngca_prepFiles.py \
    ngca_submitClusters.py \
    ngca_downloadAuspos.py \
    python3
do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        fail_controller \
            STARTUP_VALIDATION \
            "Required command not found: $required_command" \
            69
    fi
done


#
# Start the controller
#

printf '\n'
printf '%s\n' \
    '================================================================'
printf 'STARTING NGCA CONTROLLER\n'
printf 'Run ID: %s\n' "$RUN_ID"
printf 'Jurisdictions: %s\n' "${JURIS_LIST[*]^^}"
printf 'AUSPOS email: %s\n' "$EMAIL"
printf 'Submission delay: %s seconds\n' \
    "$AUSPOS_SUBMISSION_DELAY_SECONDS"
printf '%s\n' \
    '================================================================'

#
# Download and organise input data
#
# ngca_get.py is treated as a controller-level dependency because it is
# called once for all selected jurisdictions.
#

GET_OUTPUT_FILE="$(mktemp)"

cleanup_temporary_files() {
    rm -f "${GET_OUTPUT_FILE:-}"
}

trap cleanup_temporary_files EXIT

log INFO NONE DOWNLOAD \
    "Downloading and organising NGCA files"

if PYTHONUNBUFFERED=1 \
    ngca_get.py \
    -j "${JURIS_LIST[@]}" \
    > "$GET_OUTPUT_FILE" 2>&1
then
    :
else
    get_exit=$?

    cat "$GET_OUTPUT_FILE"

    fail_controller \
        DOWNLOAD \
        "ngca_get.py failed with exit code $get_exit" \
        "$get_exit"
fi

# The current ngca_get.py interface is expected to print the archive
# on the first output line.
ARCHIVE="$(head -n 1 "$GET_OUTPUT_FILE" | tr -d '\r')"

# Print all remaining downloader output.
tail -n +2 "$GET_OUTPUT_FILE" || true

if [[ ! "$ARCHIVE" =~ ^[0-9]{8}$ ]]; then
    fail_controller \
        DOWNLOAD \
        "Invalid archive returned by ngca_get.py: $ARCHIVE" \
        30
fi

log INFO NONE DOWNLOAD "NGCA download completed archive=$ARCHIVE"

NOTES_FILE="$SCRIPT_DIR/auto_notes_${ARCHIVE}.txt"

# Truncate the existing notes file or create a new one.
: > "$NOTES_FILE"

#
# Process every jurisdiction
#

for JURIS in "${JURIS_LIST[@]}"; do
    printf '\n'
    printf '%s\n' \
        '----------------------------------------------------------------'
    printf 'PROCESSING %s\n' "${JURIS^^}"
    printf '%s\n' \
        '----------------------------------------------------------------'

    if process_jurisdiction "$JURIS"; then
        printf '\n%s completed successfully.\n' "${JURIS^^}"
    else
        printf '\n%s failed at stage %s. Continuing with remaining jurisdictions.\n' \
            "${JURIS^^}" \
            "${JURIS_FAILED_STAGE[$JURIS]:-UNKNOWN}" \
            >&2
    fi
done

DOWNLOADER_EXIT=0

if run_auspos_downloader; then
    DOWNLOADER_EXIT=0
else
    DOWNLOADER_EXIT=$?
fi

#
# Final summary
#

print_summary

FAILED_COUNT=0

for JURIS in "${JURIS_LIST[@]}"; do
    if [[ "${JURIS_STATUS[$JURIS]:-UNKNOWN}" != "SUCCEEDED" ]]; then
        ((FAILED_COUNT += 1))
    fi
done

if ((FAILED_COUNT > 0)); then
    printf '\nERROR: %d jurisdiction(s) failed.\n' \
        "$FAILED_COUNT" >&2
    printf 'Downloader exit code: %d\n' \
        "$DOWNLOADER_EXIT" >&2
    
    printf 'See the CloudWatch logs for details.\n' >&2

    # AWS Batch will mark this controller job FAILED.
    # Configure the controller job with one attempt initially so Batch
    # does not automatically repeat partial AUSPOS submissions.
    exit 2
fi

printf '\nAll selected jurisdictions completed successfully.\n'
exit 0