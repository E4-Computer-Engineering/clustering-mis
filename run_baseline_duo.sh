#!/usr/bin/env bash

ml OpenMPI
ml HyperQueue

set -x
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for i in {1..5}; do
    echo "Iteration $i: Submitting jobs"

    # Submit two jobs and capture their job IDs
    # jobid1=$(sbatch --parsable --wrap="${SCRIPT_DIR}/run_baseline_single.sh perf/baseline_sleep_duo${i}_1")
    # jobid2=$(sbatch --parsable --wrap="${SCRIPT_DIR}/run_baseline_single.sh perf/baseline_sleep_duo${i}_2")
    jobid1=$(sbatch --parsable "${SCRIPT_DIR}/run_baseline_single.sh" "perf/baseline_sleep_duo${i}_1")
    jobid2=$(sbatch --parsable "${SCRIPT_DIR}/run_baseline_single.sh" "perf/baseline_sleep_duo${i}_2")
    

    echo "Submitted jobs: $jobid1 and $jobid2"

    # Function to wait for job completion
    wait_for_jobs() {
        local jobids=("$@")
        while true; do
            # Check if any job is still in the queue
            running=false
            for jobid in "${jobids[@]}"; do
                if squeue -j "$jobid" > /dev/null 2>&1 && squeue -j "$jobid" | grep -q "$jobid"; then
                    running=true
                    break
                fi
            done
            if ! $running; then
                break
            fi
            sleep 5
        done
    }

    # Wait for both jobs to finish
    wait_for_jobs "$jobid1" "$jobid2"

    echo "Jobs $jobid1 and $jobid2 finished."
done

echo "All iterations complete."

