#!/bin/bash
# Watch the last two epitopes (f4, f5) to completion.
#
# Note: `sapia run` is detached -- the tasks live on Modal, not in this script.
# The earlier campaign driver died mid-E4 and its 25 tasks carried on regardless,
# which is why f4 kept filling with no local process attached. So this only
# observes; killing it cannot stop the run.
set -u
cd /home/lpignatti/ant-comp
R=outputs/20261007_155900_tnf_1000
for i in $(seq 1 180); do
  OUT=$(NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd \
    "for L in f4 f5; do echo \"\$L \$(ls $R/table0_\$L/laproteina/laproteina_logs/*.exit 2>/dev/null | wc -l) \$(find $R/table0_\$L -name '*.pdb' 2>/dev/null | wc -l)\"; done" 2>/dev/null | tail -2)
  echo "[$(date +%H:%M)] $(echo "$OUT" | tr '\n' ' ')"
  D4=$(echo "$OUT" | awk '/^f4/{print $2}'); D5=$(echo "$OUT" | awk '/^f5/{print $2}')
  if [ "${D4:-0}" = "25" ] && [ "${D5:-0}" = "25" ]; then
    echo "=== BOTH COMPLETE $(date +%H:%M:%S) ==="
    NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd \
      "for d in $R/table0_*/; do echo \"  \$(basename \$d): \$(find \$d -name '*.pdb' 2>/dev/null | wc -l) pdbs\"; done" 2>/dev/null | tail -10
    exit 0
  fi
  sleep 120
done
