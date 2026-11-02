# Sequencer run automation using CopyComplete.txt

The sequencer writes a file called *CopyComplete.txt* to the root of each run directory when all raw data has been transferred from the instrument to its output location.
This signal file enables downstream automation such as uploading data to cloud storage or launching an analysis pipeline without manual intervention.

## Overview

A watcher script polls the sequencer output directory on a configurable interval.
It tracks the state of each discovered run directory using a hidden state directory (`.run_state/`) that contains per-run signal files.
Each run progresses through the following states:

| Signal file | Meaning |
|---|---|
| `<run_id>.found` | Run directory detected; transfer is in progress |
| `<run_id>.transfer_complete` | *CopyComplete.txt* was detected; data is ready |
| `<run_id>.timed_out` | Transfer timeout elapsed; `transfer_timeout` was called |
| `<run_id>.submitted` | Analysis has been submitted |

The state directory persists across restarts of the watcher, so runs that were already submitted are not reprocessed.

## Automation use cases

### Upload to cloud or data transfer

Once `<run_id>.transfer_complete` is signalled, the `launch_analysis` hook can be replaced or supplemented with a transfer command.
Common approaches include:

- `aws s3 sync <run_dir> s3://<bucket>/<run_id>/` for AWS S3.
- `gsutil -m rsync -r <run_dir> gs://<bucket>/<run_id>/` for Google Cloud Storage.
- `rclone copy <run_dir> remote:<bucket>/<run_id>` for any rclone-supported endpoint.

If both a transfer and an analysis are needed, run the transfer inside `launch_analysis` and wait for it to complete before submitting the pipeline.

### Launching an analysis pipeline

The `launch_analysis` function receives the absolute path to the completed run directory.
Use it to submit a job to a scheduler (SLURM, LSF) or invoke a pipeline directly.
Example for XOOS germline pipeline:

```bash
launch_analysis() {
    local run_dir="$1"
    nextflow run xoosnf/main.nf \
        --input "${run_dir}" \
        --outdir "${run_dir}/analysis" \
        -profile slurm
}
```

## Example watcher script

```bash
#!/usr/bin/env bash
set -eu -o pipefail

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") <watch-dir> [poll-interval] [transfer-timeout]" >&2
    echo "" >&2
    echo "  watch-dir          Path to the directory containing AXELIOS 1 run folders (required)" >&2
    echo "  poll-interval      Polling interval in seconds (default: 300)" >&2
    echo "  transfer-timeout   Seconds to wait for CopyComplete.txt before timing out (default: 43200)" >&2
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
WATCH_DIR="$1"
POLL_INTERVAL="${2:-300}"       # seconds; default 5 minutes
TRANSFER_TIMEOUT="${3:-43200}"  # seconds; default 12 hours
STATE_DIR="${WATCH_DIR}/.run_state"

# ---------------------------------------------------------------------------
# Placeholder functions — implement these for your environment
# ---------------------------------------------------------------------------

# Called when a run directory has been found but CopyComplete.txt has not
# appeared within a reasonable time. Receives the run directory as $1.
# Use this to alert, log, or escalate a stalled transfer.
transfer_timeout() {
    local run_dir="$1"
    echo "[WARN] Transfer timeout for run: ${run_dir}" >&2
    # TODO: implement alerting (e.g., send email, post to Slack, page on-call)
}

# Called once CopyComplete.txt is detected in a run directory. Receives the
# run directory as $1. Launch whatever downstream analysis is required.
launch_analysis() {
    local run_dir="$1"
    echo "[INFO] Launching analysis for run: ${run_dir}"
    # TODO: implement analysis submission, for example:
    #   nextflow run xoosnf/main.nf --input "${run_dir}" -profile slurm
    #   aws s3 sync "${run_dir}" s3://my-bucket/"$(basename "${run_dir}")"
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

state_file() {
    local run_id="$1"
    local state="$2"
    echo "${STATE_DIR}/${run_id}.${state}"
}

has_state() {
    local run_id="$1"
    local state="$2"
    [[ -f "$(state_file "${run_id}" "${state}")" ]]
}

set_state() {
    local run_id="$1"
    local state="$2"
    touch "$(state_file "${run_id}" "${state}")"
    echo "[INFO] $(date -Iseconds) run=${run_id} state=${state}"
}

# Returns the number of seconds elapsed since the .found signal file was written.
found_age_seconds() {
    local run_id="$1"
    local found_file
    found_file="$(state_file "${run_id}" "found")"
    echo $(( $(date +%s) - $(date +%s -r "${found_file}") ))
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

mkdir -p "${STATE_DIR}"
echo "[INFO] Watching ${WATCH_DIR} every ${POLL_INTERVAL}s"

while true; do
    for run_dir in "${WATCH_DIR}"/*/; do
        # Skip if glob matched nothing
        [[ -d "${run_dir}" ]] || continue

        run_id="$(basename "${run_dir}")"

        # Skip the state directory itself
        [[ "${run_id}" == ".run_state" ]] && continue

        # Step 1: Mark newly discovered runs as found
        if ! has_state "${run_id}" "found"; then
            set_state "${run_id}" "found"
        fi

        # Step 2: Check whether data transfer is complete
        if ! has_state "${run_id}" "transfer_complete"; then
            if [[ -f "${run_dir}/CopyComplete.txt" ]]; then
                set_state "${run_id}" "transfer_complete"
            elif (( $(found_age_seconds "${run_id}") > TRANSFER_TIMEOUT )); then
                if ! has_state "${run_id}" "timed_out"; then
                    set_state "${run_id}" "timed_out"
                    transfer_timeout "${run_dir}"
                fi
                continue
            else
                continue
            fi
        fi

        # Step 3: Submit analysis for runs that are complete but not yet submitted
        if ! has_state "${run_id}" "submitted"; then
            launch_analysis "${run_dir}"
            set_state "${run_id}" "submitted"
        fi
    done

    sleep "${POLL_INTERVAL}"
done
```

### Running the script

```bash
# Watch a specific directory with default poll interval (5 min) and timeout (12 hr)
bash watch_sequencer.sh /data/sequencer_output

# Override poll interval (60 s) and transfer timeout (24 hr)
bash watch_sequencer.sh /mnt/instrument_output 60 86400

# Run as a systemd service or in the background with logging
nohup bash watch_sequencer.sh /data/sequencer_output >> /var/log/sequencer_watcher.log 2>&1 &
```

### State directory layout

After several runs the state directory looks like this:

```text
/data/sequencer_output/
├── .run_state/
│   ├── 20260601_run001.found
│   ├── 20260601_run001.transfer_complete
│   ├── 20260601_run001.submitted
│   ├── 20260602_run002.found
│   └── 20260602_run002.transfer_complete   ← analysis not yet submitted
├── 20260601_run001/
│   ├── CopyComplete.txt
│   └── ...
└── 20260602_run002/
    ├── CopyComplete.txt
    └── ...
```

## Known issues

### CopyComplete.txt may not be created in rare cases

In rare circumstances *CopyComplete.txt* may never appear even though the run directory is otherwise complete (for example, after a software crash or power loss).
The transfer timeout implemented in the script will detect this and call `transfer_timeout` so operators can intervene.

### No indication of sequencing error in the signal file

*CopyComplete.txt* confirms that the instrument finished writing data — it does not indicate whether the run produced sufficient or valid data.
There is no clear automated workaround; always check alignment and yield metrics after analysis to confirm sequencing output is inline with expectations.
