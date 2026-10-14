#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=seqcluster

# Cluster the WHOLE run set in this one task (seqcluster submits exactly one task per
# run: clustering is a set-level operation). The manifest line points at a members
# file (name<TAB>sequence per design) and carries the run token and the four MMseqs2
# parameters. The worker writes seqcluster_result/assignments.tsv, one row per design,
# for collect_seqcluster.py. CPU-only. The worker is stdlib-only Python and shells out
# to `mmseqs` ($MMSEQS, default `mmseqs` on PATH).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_SEQCLUSTER

PY=${PIPELINE_PYTHON:-python}
MMSEQS_BIN=${MMSEQS:-mmseqs}

# NOTE: never name a manifest field after a bash special variable (GROUPS, UID,
# PIPESTATUS, ...) — the assignment fails with rc=1 and, under `set -e`, the task
# dies with EMPTY .out and .err.
MEMBERS_FILE=$(echo "$SAPIA_LINE" | cut -f1)
RUN_TOKEN=$(echo "$SAPIA_LINE" | cut -f2)
MIN_SEQ_ID=$(echo "$SAPIA_LINE" | cut -f3)
COVERAGE=$(echo "$SAPIA_LINE" | cut -f4)
COV_MODE=$(echo "$SAPIA_LINE" | cut -f5)
CLUSTER_MODE=$(echo "$SAPIA_LINE" | cut -f6)
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f7)

echo "[$(date +%T)] task $SAPIA_TASK_ID: seqcluster on $(wc -l <"$MEMBERS_FILE") designs" \
     "(min_seq_id=$MIN_SEQ_ID c=$COVERAGE cov_mode=$COV_MODE cluster_mode=$CLUSTER_MODE)"

# Clear the previous result BEFORE the worker runs, so a task killed mid-run (OOM,
# timeout, pre-emption) leaves NO result rather than an earlier run's. The run token
# the collector checks is the real guard; this keeps the directory honest too.
rm -f "$OUT_DIR/seqcluster_result/assignments.tsv"

# No per-design fallback rows: if the worker dies without writing assignments.tsv,
# the collector marks every member 'error:' from its own copy of the members list.
# The worker's return code is captured EXPLICITLY and re-raised, so the task's exit
# status means what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/seqcluster_worker.py" \
        --members "$MEMBERS_FILE" \
        --token "$RUN_TOKEN" \
        --min-seq-id "$MIN_SEQ_ID" \
        --coverage "$COVERAGE" \
        --cov-mode "$COV_MODE" \
        --cluster-mode "$CLUSTER_MODE" \
        --input-column "$INPUT_COLUMN" \
        --out-dir "$OUT_DIR" \
        --mmseqs "$MMSEQS_BIN" || WORKER_RC=$?

exit "$WORKER_RC"
