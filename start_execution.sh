#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --output=output_%j.txt   
#SBATCH --time=01:30:00         
#SBATCH --nodes=1            

dmr_wrapper mpirun -n 1 build/bin/clustering data/input/cluster_points_article.csv work_dir

