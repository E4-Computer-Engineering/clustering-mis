#!/usr/bin/env bash

set -x

ml OpenMPI
source .venv/bin/activate

for i in {1..5};
do
    mkdir -p "perf/workflow${i}"
    { time streamflow run --name "smart-hpc-qc${i}" workflow/streamflow-e4.yml; } > "perf/workflow${i}/log.txt" 2>&1
    streamflow report --format csv --file workflow/streamflow-e4.yml "smart-hpc-qc${i}"
    mv report.csv "perf/workflow${i}/"
done
