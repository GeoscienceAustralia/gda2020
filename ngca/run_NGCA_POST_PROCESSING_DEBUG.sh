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
ARCHIVE=""
RUN_ID=""

CPUS="${NGCA_POSTPROCESS_CPUS:-4}"

# Optional. When empty, existing local files are used.
NGCA_RUN_BUCKET="${NGCA_RUN_BUCKET:-}"
NGCA_RESULT_DEST="${NGCA_RESULT_DEST:-}"

declare -A JURIS_STATUS
declare -A JURIS_FAILED_STAGE
declare -A JURIS_MESSAGE

usage() {
    cat <<EOF
Usage:
  $(basename "$0") \
      --juris JURIS [JURIS ...] \
      --archive YYYYMMDD \
      --run-id RUN_ID \
      [--cpu CPUS]

Options:
  -j, --juris JURIS...   Jurisdictions to process
  -a, --archive VALUE    NGCA archive date, for example 20261105
  -r, --run-id VALUE     NGCA run identifier
  -c, --cpu VALUE        CPUs for minimal adjustment
                         Default: NGCA_POSTPROCESS_CPUS or 4
  -h, --help             Show this help

Environment variables:
  NGCA_RUN_BUCKET
      Optional S3 bucket containing the post-processing input zip.
      If unset, existing local files are used.

  NGCA_POSTPROCESS_CPUS
      Default CPU count when --cpu is not supplied.

Expected S3 object:
  s3://<bucket>/NGCA/batch_jobs/<juris>/<run-id>/NGCA-<JURIS>-<run-id>.zip

Expected zip contents:
  <archive>/
  rinexantls/
  sinexFiles/
EOF
}

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

add_jurisdiction() {
    local requested="${1,,}"
    local existing

    if ! is_valid_jurisdiction "$requested"; then
        printf 'ERROR: Unknown jurisdiction: %s\n' \
            "$requested" >&2
        exit 64
    fi

    for existing in "${JURIS_LIST[@]}"; do
        if [[ "$existing" == "$requested" ]]; then
            return 0
        fi
    done

    JURIS_LIST+=("$requested")
}

mark_jurisdiction_running() {
    local jurisdiction="$1"

    JURIS_STATUS["$jurisdiction"]="RUNNING"
    JURIS_FAILED_STAGE["$jurisdiction"]=""
    JURIS_MESSAGE["$jurisdiction"]=""

    log INFO "$jurisdiction" START \
        "Post-processing started"
}

mark_jurisdiction_failed() {
    local jurisdiction="$1"
    local stage="$2"
    local message="$3"

    JURIS_STATUS["$jurisdiction"]="FAILED"
    JURIS_FAILED_STAGE["$jurisdiction"]="$stage"
    JURIS_MESSAGE["$jurisdiction"]="$message"

    log ERROR "$jurisdiction" "$stage" "$message"
}

mark_jurisdiction_complete() {
    local jurisdiction="$1"

    JURIS_STATUS["$jurisdiction"]="SUCCEEDED"
    JURIS_FAILED_STAGE["$jurisdiction"]=""
    JURIS_MESSAGE["$jurisdiction"]=""

    log INFO "$jurisdiction" COMPLETE \
        "NGCA post-processing completed successfully"
}

# Return s3 key to zip file
postprocess_object_key() {
    local jurisdiction="$1"

    printf 'NGCA/batch_jobs/%s/%s/NGCA-%s-%s.zip' \
        "$jurisdiction" \
        "$RUN_ID" \
        "${jurisdiction^^}" \
        "$RUN_ID"
}

# Download zip file from s3 bucket if bucket is given
restore_postprocessing_input() {
    local jurisdiction="$1"
    local jurisdiction_directory
    local object_key
    local zip_file
    local temporary_directory

    # If no bucket provided then use local files
    if [[ -z "$NGCA_RUN_BUCKET" ]]; then
        log INFO "$jurisdiction" INPUT_RESTORE \
            "No NGCA_RUN_BUCKET configured; using local files"

        return 0
    fi

    jurisdiction_directory="$SCRIPT_DIR/$jurisdiction"
    object_key="$(postprocess_object_key "$jurisdiction")"

    temporary_directory="$(mktemp -d)"
    zip_file="$temporary_directory/NGCA-${jurisdiction^^}-${RUN_ID}.zip"

    log INFO "$jurisdiction" INPUT_RESTORE \
        "Downloading s3://${NGCA_RUN_BUCKET}/${object_key}"

    # Try to download to temporary directory
    if ! aws s3 cp \
        "s3://${NGCA_RUN_BUCKET}/${object_key}" \
        "$zip_file"
    then
        rm -rf "$temporary_directory"

        log ERROR "$jurisdiction" INPUT_RESTORE \
            "Could not download post-processing input from S3"

        return 1
    fi

    if [[ ! -s "$zip_file" ]]; then
        rm -rf "$temporary_directory"

        log ERROR "$jurisdiction" INPUT_RESTORE \
            "Downloaded post-processing zip is empty"

        return 1
    fi

    if ! unzip -t "$zip_file" >/dev/null; then
        rm -rf "$temporary_directory"

        log ERROR "$jurisdiction" INPUT_RESTORE \
            "Downloaded post-processing zip is invalid"

        return 1
    fi

    log INFO "$jurisdiction" INPUT_RESTORE \
        "Extracting post-processing input into $jurisdiction_directory"

    # Try to unzip file
    if ! unzip -o \
        "$zip_file" \
        -d "$jurisdiction_directory"
    then
        rm -rf "$temporary_directory"

        log ERROR "$jurisdiction" INPUT_RESTORE \
            "Could not extract post-processing input"

        return 1
    fi

    rm -rf "$temporary_directory"

    log INFO "$jurisdiction" INPUT_RESTORE \
        "Post-processing input restored successfully"
}

validate_input() {
    local jurisdiction="$1"

    local jurisdiction_directory
    local archive_directory
    local antenna_directory
    local sinex_directory

    jurisdiction_directory="$SCRIPT_DIR/$jurisdiction"
    archive_directory="$jurisdiction_directory/$ARCHIVE"
    antenna_directory="$jurisdiction_directory/rinexantls"
    sinex_directory="$jurisdiction_directory/sinexFiles"

    if [[ ! -d "$archive_directory" ]]; then
        log ERROR "$jurisdiction" INPUT_VALIDATION \
            "Archive directory not found: $archive_directory"
        return 1
    fi

    if [[ ! -d "$antenna_directory" ]]; then
        log ERROR "$jurisdiction" INPUT_VALIDATION \
            "RINEX antenna directory not found: $antenna_directory"
        return 1
    fi

    if [[ ! -d "$sinex_directory" ]]; then
        log ERROR "$jurisdiction" INPUT_VALIDATION \
            "SINEX directory not found: $sinex_directory"
        return 1
    fi

    if ! find "$sinex_directory" \
        -maxdepth 1 \
        -type f \
        -iname '*.snx' \
        -print \
        -quit |
        grep -q .
    then
        log ERROR "$jurisdiction" INPUT_VALIDATION \
            "No SINEX files found in $sinex_directory"
        return 1
    fi

    log INFO "$jurisdiction" INPUT_VALIDATION \
        "Post-processing input validation passed"
}

run_stage() {
    local jurisdiction="$1"
    local stage="$2"
    local description="$3"

    shift 3

    log INFO "$jurisdiction" "$stage" "$description"

    if "$@"; then
        log INFO "$jurisdiction" "$stage" \
            "$description completed successfully"

        return 0
    fi

    local exit_code=$?

    mark_jurisdiction_failed \
        "$jurisdiction" \
        "$stage" \
        "$description failed with exit code $exit_code"

    return "$exit_code"
}

process_jurisdiction() {
    local jurisdiction="$1"
    local archive_directory

    archive_directory="$SCRIPT_DIR/$jurisdiction/$ARCHIVE"

    mark_jurisdiction_running "$jurisdiction"

    if ! restore_postprocessing_input "$jurisdiction"; then
        mark_jurisdiction_failed \
            "$jurisdiction" \
            INPUT_RESTORE \
            "Could not restore post-processing input"

        return 1
    fi

    if ! validate_input "$jurisdiction"; then
        mark_jurisdiction_failed \
            "$jurisdiction" \
            INPUT_VALIDATION \
            "Post-processing input validation failed"

        return 1
    fi

    cd "$archive_directory"

    if ! run_stage \
        "$jurisdiction" \
        TRANSFORM_SINEX \
        "Transforming SINEX from ITRF2020 to ITRF2014" \
        ngca_transformSinex.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        EXCISE_STATIONS \
        "Removing stations outside GDA2020" \
        ngca_exciseStations.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        SELECT_REFERENCE \
        "Selecting APREF constraint stations" \
        ngca_selectRef.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        CREATE_BASELINES \
        "Creating baselines and DynaML files" \
        ngca_createBLs.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        UPDATE_TRANSLATION_TABLES \
        "Updating translation tables" \
        ngca_updateTTables.py
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        CHECK_TRANSLATION_TABLES \
        "Checking translation tables" \
        ngca_checkTTables.py
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        TRANSLATE_STATIONS \
        "Updating station names using translation tables" \
        ngca_stationTrans.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        MINIMAL_ADJUSTMENT \
        "Running minimally constrained cluster adjustments" \
        ngca_minAdjust.py \
            --juris "$jurisdiction" \
            --cpu "$CPUS"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        UPDATE_VSCALE \
        "Updating Vscale values from sigma0.dat" \
        ngca_addSigma0.py \
            --juris "$jurisdiction"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        COMBINE_CLUSTERS \
        "Combining adjusted clusters" \
        ngca_combine.py \
            --juris "$jurisdiction" \
            --archive "$ARCHIVE"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        FULL_ADJUSTMENT \
        "Running fully constrained adjustment" \
        ngca_fullAdjust.py \
            --juris "$jurisdiction" \
            --archive "$ARCHIVE"
    then
        return 1
    fi

    if ! run_stage \
        "$jurisdiction" \
        ORGANISE_RESULTS \
        "Organising NGCA results" \
        ngca_results.py \
            --juris "$jurisdiction" \
            --archive "$ARCHIVE" \
            --destination "$NGCA_RESULT_DEST"
    then
        return 1
    fi

    mark_jurisdiction_complete "$jurisdiction"
    return 0
}

print_summary() {
    local jurisdiction
    local status
    local passed=0
    local failed=0

    printf '\n'
    printf '%s\n' \
        '================================================================'
    printf 'NGCA POST-PROCESSING SUMMARY\n'
    printf 'Run ID:  %s\n' "$RUN_ID"
    printf 'Archive: %s\n' "$ARCHIVE"
    printf 'CPUs:    %s\n' "$CPUS"
    printf '%s\n' \
        '================================================================'

    for jurisdiction in "${JURIS_LIST[@]}"; do
        status="${JURIS_STATUS[$jurisdiction]:-UNKNOWN}"

        case "$status" in
            SUCCEEDED)
                ((passed += 1))

                printf '%-4s  %-10s\n' \
                    "${jurisdiction^^}" \
                    "SUCCEEDED"
                ;;

            FAILED)
                ((failed += 1))

                printf '%-4s  %-10s  STAGE=%-28s  %s\n' \
                    "${jurisdiction^^}" \
                    "FAILED" \
                    "${JURIS_FAILED_STAGE[$jurisdiction]:-UNKNOWN}" \
                    "${JURIS_MESSAGE[$jurisdiction]:-No message}"
                ;;

            *)
                ((failed += 1))

                printf '%-4s  %-10s  %s\n' \
                    "${jurisdiction^^}" \
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
# Parse arguments
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

        -a|--archive)
            if (($# < 2)); then
                printf 'ERROR: %s requires a value\n' "$1" >&2
                exit 64
            fi

            ARCHIVE="$2"
            shift 2
            ;;

        -r|--run-id)
            if (($# < 2)); then
                printf 'ERROR: %s requires a value\n' "$1" >&2
                exit 64
            fi

            RUN_ID="$2"
            shift 2
            ;;

        -c|--cpu)
            if (($# < 2)) || [[ ! "$2" =~ ^[1-9][0-9]*$ ]]; then
                printf 'ERROR: --cpu requires a positive integer\n' \
                    >&2
                exit 64
            fi

            CPUS="$2"
            shift 2
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
            printf 'ERROR: Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 64
            ;;
    esac
done

#
# Validate arguments
#

if ((${#JURIS_LIST[@]} == 0)); then
    printf 'ERROR: At least one jurisdiction is required\n' \
        >&2
    exit 64
fi

if [[ ! "$ARCHIVE" =~ ^[0-9]{8}$ ]]; then
    printf 'ERROR: --archive must use YYYYMMDD format\n' \
        >&2
    exit 64
fi

if [[ -z "$RUN_ID" ]]; then
    printf 'ERROR: --run-id is required\n' >&2
    exit 64
fi

if [[ ! "$RUN_ID" =~ ^[a-zA-Z0-9][a-zA-Z0-9-]*$ ]]; then
    printf 'ERROR: Invalid run ID: %s\n' "$RUN_ID" >&2
    printf 'Use letters, numbers and hyphens only.\n' >&2
    exit 64
fi

#
# Configure script path
#

if [[ -d "$SCRIPT_DIR/scripts" ]] &&
   [[ ":$PATH:" != *":$SCRIPT_DIR/scripts:"* ]]; then
    export PATH="$SCRIPT_DIR/scripts:$PATH"
fi

#
# Validate dependencies
#

required_commands=(
    aws
    unzip
    grep
    find
    ngca_transformSinex.py
    ngca_exciseStations.py
    ngca_selectRef.py
    ngca_createBLs.py
    ngca_updateTTables.py
    ngca_checkTTables.py
    ngca_stationTrans.py
    ngca_minAdjust.py
    ngca_addSigma0.py
    ngca_combine.py
    ngca_fullAdjust.py
    ngca_results.py
)

for required_command in "${required_commands[@]}"; do
    # AWS CLI is optional when no S3 bucket is configured.
    if [[ "$required_command" == "aws" ]] &&
       [[ -z "$NGCA_RUN_BUCKET" ]]; then
        continue
    fi

    if ! command -v "$required_command" >/dev/null 2>&1; then
        printf 'ERROR: Required command not found: %s\n' \
            "$required_command" >&2
        exit 69
    fi
done

#
# Start processing
#

printf '\n'
printf '%s\n' \
    '================================================================'
printf 'STARTING NGCA POST-PROCESSING\n'
printf 'Run ID:        %s\n' "$RUN_ID"
printf 'Archive:       %s\n' "$ARCHIVE"
printf 'Jurisdictions: %s\n' "${JURIS_LIST[*]^^}"
printf 'CPUs:          %s\n' "$CPUS"

if [[ -n "$NGCA_RUN_BUCKET" ]]; then
    printf 'Input source:   s3://%s\n' "$NGCA_RUN_BUCKET"
else
    printf 'Input source:   local files\n'
fi

printf '%s\n' \
    '================================================================'

for jurisdiction in "${JURIS_LIST[@]}"; do
    printf '\n'
    printf '%s\n' \
        '----------------------------------------------------------------'
    printf 'POST-PROCESSING %s\n' "${jurisdiction^^}"
    printf '%s\n' \
        '----------------------------------------------------------------'

    if process_jurisdiction "$jurisdiction"; then
        printf '\n%s post-processing completed successfully.\n' \
            "${jurisdiction^^}"
    else
        printf '\n%s post-processing failed at stage %s.\n' \
            "${jurisdiction^^}" \
            "${JURIS_FAILED_STAGE[$jurisdiction]:-UNKNOWN}" \
            >&2
    fi
done

print_summary

failed_count=0

for jurisdiction in "${JURIS_LIST[@]}"; do
    if [[ "${JURIS_STATUS[$jurisdiction]:-UNKNOWN}" != "SUCCEEDED" ]]; then
        ((failed_count += 1))
    fi
done

if ((failed_count > 0)); then
    printf '\nERROR: %d jurisdiction(s) failed post-processing.\n' \
        "$failed_count" >&2

    exit 2
fi

printf '\nAll selected jurisdictions completed post-processing successfully.\n'
exit 0