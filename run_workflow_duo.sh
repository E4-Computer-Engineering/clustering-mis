#!/usr/bin/env bash

ml OpenMPI
source .venv/bin/activate

set -x

if [ ! -f "workflow/streamflow-e4-bis.yml" ]; then
  echo "Error: file 'workflow/streamflow-e4-bis.yml' does not exist."
  echo "You can copy it from the regular workflow/streamflow-e4.yml file, but"
  echo "make sure to change the database file (edit the 'connection' value)"
  exit 1
fi

# NOTE: the two loops have been actually run in two different moments. The SimulatedAnnealing executable was different (no sleep vs 2m sleep).

# SimulatedAnnealing finishes after a fraction of a second
for i in {1..5};
do
    mkdir -p "perf/workflow_duo${i}_1"
    mkdir -p "perf/workflow_duo${i}_2"

    { time streamflow run --name "smart-hpc-qc_duo${i}_1" workflow/streamflow-e4.yml; } > "perf/workflow_duo${i}_1/log.txt" 2>&1 &
    { time streamflow run --name "smart-hpc-qc_duo${i}_2" workflow/streamflow-e4-bis.yml; } > "perf/workflow_duo${i}_2/log.txt" 2>&1 &
    
    wait

    streamflow report --format csv --file workflow/streamflow-e4.yml "smart-hpc-qc_duo${i}_1"
    mv report.csv "perf/workflow_duo${i}_1/"

    streamflow report --format csv --file workflow/streamflow-e4-bis.yml "smart-hpc-qc_duo${i}_2"
    mv report.csv "perf/workflow_duo${i}_2/"
    
done

# SimulatedAnnealing sleeps for 2 minutes before giving its result
for i in {1..5};
do
    mkdir -p "perf/workflow_sleep_duo${i}_1"
    mkdir -p "perf/workflow_sleep_duo${i}_2"

    { time streamflow run --name "smart-hpc-qc_sleep_duo${i}_1" workflow/streamflow-e4.yml; } > "perf/workflow_sleep_duo${i}_1/log.txt" 2>&1 &
    { time streamflow run --name "smart-hpc-qc_sleep_duo${i}_2" workflow/streamflow-e4-bis.yml; } > "perf/workflow_sleep_duo${i}_2/log.txt" 2>&1 &
    
    wait

    streamflow report --format csv --file workflow/streamflow-e4.yml "smart-hpc-qc_sleep_duo${i}_1"
    mv report.csv "perf/workflow_sleep_duo${i}_1/"

    streamflow report --format csv --file workflow/streamflow-e4-bis.yml "smart-hpc-qc_sleep_duo${i}_2"
    mv report.csv "perf/workflow_sleep_duo${i}_2/"
    
done
