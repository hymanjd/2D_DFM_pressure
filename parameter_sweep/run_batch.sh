#!/usr/bin/env bash
# run_sweep.sh — parallel DFM sweep, N jobs at a time
# Usage: bash run_sweep.sh [njobs]

NJOBS=${1:-8}
N=$(python3.11 -c "import pandas as pd; print(len(pd.read_csv('lhs_params.csv')))")
START=$(date +%s)

echo "=================================================="
echo " DFM sweep: $N samples  |  $NJOBS parallel slots"
echo " Start: $(date)"
echo "=================================================="

declare -A JOB_START   # pid -> start epoch
declare -A JOB_IDX     # pid -> sample index
running=0
finished=0
failed=0

finish_one() {
    local pid=$1
    wait "$pid"
    local rc=$?
    local idx=${JOB_IDX[$pid]}
    local elapsed=$(( $(date +%s) - ${JOB_START[$pid]} ))
    local mins=$(( elapsed / 60 ))
    local secs=$(( elapsed % 60 ))
    if [[ $rc -eq 0 ]]; then
        (( finished++ ))
        printf "  [%04d] done   %dm%02ds   (%d/%d finished, %d failed)\n" \
            "$idx" "$mins" "$secs" "$finished" "$N" "$failed"
    else
        (( failed++ ))
        printf "  [%04d] FAILED %dm%02ds   (rc=%d)  -- see x%04d.out\n" \
            "$idx" "$mins" "$secs" "$rc" "$idx"
    fi
    unset JOB_START[$pid]
    unset JOB_IDX[$pid]
}

for i in $(seq 1 $N); do
    outfile=$(printf "x%04d.out" "$i")
    python3.11 dfm_driver.py "$i" > "$outfile" 2>&1 &
    pid=$!
    JOB_START[$pid]=$(date +%s)
    JOB_IDX[$pid]=$i
    (( running++ ))
    printf "  [%04d] launched  (pid %d)  [%d running]\n" "$i" "$pid" "$running"

    if (( running >= NJOBS )); then
        # wait for any one child to finish
        wait -n 2>/dev/null
        exited_pid=$!
        # bash 5.1+: wait -n -p exited_pid; older bash needs the loop below
        if [[ -z "${JOB_IDX[$exited_pid]+x}" ]]; then
            # fallback: scan for any pid that has exited
            for pid in "${!JOB_IDX[@]}"; do
                if ! kill -0 "$pid" 2>/dev/null; then
                    exited_pid=$pid; break
                fi
            done
        fi
        finish_one "$exited_pid"
        (( running-- ))
    fi
done

# drain remaining jobs
for pid in "${!JOB_IDX[@]}"; do
    finish_one "$pid"
done
wait

ELAPSED=$(( $(date +%s) - START ))
echo "=================================================="
echo " Sweep complete: $finished ok, $failed failed"
printf " Wall time: %dh %dm %02ds\n" \
    $(( ELAPSED/3600 )) $(( (ELAPSED%3600)/60 )) $(( ELAPSED%60 ))
echo " End: $(date)"
echo "=================================================="


