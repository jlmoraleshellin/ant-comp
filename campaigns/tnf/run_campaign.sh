#!/bin/bash
# TNF-alpha binder campaign: 5 epitopes x 200 designs with La-Proteina Complexa.
#
# 25 jobs x 8 designs rather than 1 job x 200. A job's dataset is 8 samples, so
# the DataLoader batch cannot exceed 8 whatever generation.dataloader.batch_size
# resolves to in the container -- which matters because an explicit
# ++generation.dataloader.batch_size=8 was accepted on the command line and
# still came back as 16 inside the task, OOM-ing a 40 GiB A100 on a ~520-residue
# complex. Constraining the dataset size is the one lever that cannot be ignored.
#
# -C 1 keeps this to a single GPU at a time.
set -u
cd /home/lpignatti/ant-comp
S=/tmp/claude-1000/-home-lpignatti/9763510a-1717-427b-a955-6e5a469c2950/scratchpad
R=outputs/20261007_155900_tnf_1000
SUM=""

for SPEC in "1:f1" "2:f2" "3:f3" "4:f4" "5:f5"; do
  E=${SPEC%%:*}; L=${SPEC##*:}
  echo "######## EPITOPE E$E (label $L) ######## $(date +%H:%M:%S)"

  NO_COLOR=1 .venv/bin/sapia modal-shell --cmd \
    "sapia run laproteina $R --task-name TNF_E$E --num-jobs 25 --nsamples 200 --design-prefix $L --table-label $L --set generation.dataloader.batch_size=8 -e modal -C 1" \
    > "$S/run_$L.log" 2>&1

  if grep -q ImageBuildError "$S/run_$L.log"; then
    echo "  E$E: BUILD FAILED"; SUM="$SUM\nE$E: BUILD FAILED"; continue
  fi
  NO_COLOR=1 grep -E "Submitting .* task" "$S/run_$L.log" | head -1

  for i in $(seq 1 150); do
    sleep 60
    D=$(NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd \
        "ls $R/table0_$L/laproteina/laproteina_logs/*.exit 2>/dev/null | wc -l" 2>/dev/null | tail -1)
    P=$(NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd \
        "find $R/table0_$L -name '*.pdb' 2>/dev/null | wc -l" 2>/dev/null | tail -1)
    D=${D:-0}; P=${P:-0}
    if [ $((i % 5)) -eq 0 ]; then echo "   [$(date +%H:%M)] tasks ${D}/25  pdbs ${P}"; fi
    # Bail out early if the first handful of tasks all failed without producing
    # anything -- no point spending the night on 25 copies of the same error.
    if [ "$D" -ge 3 ] && [ "$P" = "0" ]; then
      echo "  E$E: first $D tasks produced no structures - ABORTING epitope"
      SUM="$SUM\nE$E: ABORTED (first $D tasks produced nothing)"
      break
    fi
    if [ "$D" = "25" ]; then
      OK=$(NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd \
           "cat $R/table0_$L/laproteina/laproteina_logs/*.exit 2>/dev/null | grep -c '^0$'" 2>/dev/null | tail -1)
      echo "  E$E done after ~${i} min: ${OK:-?}/25 tasks ok, $P pdbs"
      SUM="$SUM\nE$E: ${OK:-?}/25 tasks ok, $P pdbs"
      break
    fi
  done
done

echo "================ CAMPAIGN SUMMARY $(date +%H:%M:%S) ================"
printf "%b\n" "$SUM"
TOT=$(NO_COLOR=1 timeout 280 .venv/bin/sapia modal-shell --cmd "find $R -name '*.pdb' 2>/dev/null | wc -l" 2>/dev/null | tail -1)
echo "TOTAL pdbs across campaign: ${TOT:-?}"
echo "run_dir: $R"
