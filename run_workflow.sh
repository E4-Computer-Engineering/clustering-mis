#!/usr/bin/env bash

ml OpenMPI
source .venv/bin/activate

set -x

# NOTE: the two loops have been actually run in two different moments. The SimulatedAnnealing executable was different (no sleep vs 2m sleep).

# SimulatedAnnealing finishes after a fraction of a second
for i in {1..5};
do
    mkdir -p "perf/workflow${i}"
    { time streamflow run --name "smart-hpc-qc${i}" workflow/streamflow-e4.yml; } > "perf/workflow${i}/log.txt" 2>&1
    streamflow report --format csv --file workflow/streamflow-e4.yml "smart-hpc-qc${i}"
    mv report.csv "perf/workflow${i}/"
done

# SimulatedAnnealing sleeps for 2 minutes before giving its result
for i in {1..5};
do
    mkdir -p "perf/workflow_sleep${i}"
    { time streamflow run --name "smart-hpc-qc_sleep${i}" workflow/streamflow-e4.yml; } > "perf/workflow_sleep${i}/log.txt" 2>&1
    streamflow report --format csv --file workflow/streamflow-e4.yml "smart-hpc-qc_sleep${i}"
    mv report.csv "perf/workflow_sleep${i}/"
done
