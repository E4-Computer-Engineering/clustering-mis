#!/bin/bash
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --job-name=clustering_dyn
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
##SBATCH --nodelist=awnode05
##SBATCH --wait

module load HyperQueue

home=/home/agarcia/clustering-mis

# Variables para que Hydra use SSH y pase el entorno correcto
export I_MPI_HYDRA_BOOTSTRAP=ssh
export I_MPI_HYDRA_ENV_LIST=PATH,LD_LIBRARY_PATH

# Opcional: evitar preguntas de SSH known_hosts (ajusta si tu cluster ya tiene confianza SSH configurada)
export I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS="-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"


#mpiexec.flexmpi --dynamic --excl-node ${home}/build/bin/clustering ${home}/data/input/cluster_points_article.csv ${home}/work_dir
echo "JobID: $SLURM_JOB_ID"
date
mpiexec.flexmpi --dynamic --excl-node ${home}/build/bin/clustering ${home}/data/input/cluster_points_80k.csv ${home}/work_dir
echo "JobID: $SLURM_JOB_ID"
date
