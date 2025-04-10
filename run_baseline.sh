#!/usr/bin/env bash
#SBATCH --job-name=loop_no_sleep
#SBATCH --partition=compute
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

set -x

ml OpenMPI
ml HyperQueue

for i in {1..5};
do
    mkdir -p "perf/baseline${i}"
    { time mpirun build/bin/clustering data/input/cluster_points_80k.csv "perf/baseline${i}"; } > "perf/baseline${i}/log.txt" 2>&1
done
