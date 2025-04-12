#!/usr/bin/env bash
#SBATCH --job-name=loop_no_sleep
#SBATCH --partition=compute
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

folder=$1
mkdir -p "$folder"
{ time mpirun build/bin/clustering data/input/cluster_points_80k.csv "$folder"; } > "$folder/log.txt" 2>&1
