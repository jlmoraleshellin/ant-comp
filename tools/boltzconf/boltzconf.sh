#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=boltzconf

# One batch per task; line: sub-manifest binder target contact_cutoff input_col ipsae_pae_cutoff

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_BOLTZCONF

PY=${PIPELINE_PYTHON:-python}

# Never name a field after a bash special (GROUPS, UID...): under `set -e` the task dies with empty .out.
TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
BINDER_CHAINS=$(echo "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f3)
CONTACT_CUTOFF=$(echo "$SAPIA_LINE" | cut -f4)
CONTACT_CUTOFF=${CONTACT_CUTOFF:-5.0}
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f5)
INPUT_COLUMN=${INPUT_COLUMN:-the input structure column}
# Field 6 (older manifests lack it): KEEP the 15.0 default in sync with boltzconf_worker.py.
IPSAE_PAE_CUTOFF=$(echo "$SAPIA_LINE" | cut -f6)
IPSAE_PAE_CUTOFF=${IPSAE_PAE_CUTOFF:-15.0}

echo "[$(date +%T)] task $SAPIA_TASK_ID: boltzconf on $(wc -l <"$TASK_FILE") designs" \
     "(binder=$BINDER_CHAINS target=$TARGET_CHAINS contact=$CONTACT_CUTOFF A" \
     "ipsae_pae=$IPSAE_PAE_CUTOFF A)"

# A worker crash writes a fallback error TSV per design without one; its rc is re-raised so .exit is honest.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/boltzconf_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chains "$BINDER_CHAINS" \
        --target-chains "$TARGET_CHAINS" \
        --contact-cutoff "$CONTACT_CUTOFF" \
        --ipsae-pae-cutoff "$IPSAE_PAE_CUTOFF" \
        --input-column "$INPUT_COLUMN" \
        --out-dir "$OUT_DIR" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    while IFS=$'\t' read -r NAME _SRC || [ -n "$NAME" ]; do
        [ -n "$NAME" ] || continue
        if [ ! -f "$OUT_DIR/${NAME}.tsv" ]; then
            printf 'name\tstatus\tpath\n%s\terror: worker crashed\t\n' "$NAME" \
                >"$OUT_DIR/${NAME}.tsv"
        fi
    done <"$TASK_FILE"
fi

exit "$WORKER_RC"
