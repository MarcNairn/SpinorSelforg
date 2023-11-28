#!/bin/bash
# Number of nodes to allocate
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node
#SBATCH --ntasks-per-node=1
# Number of threads per MPI instance
#SBATCH --cpus-per-task=1
# Allocate xx GB memory per node
#SBATCH --mem=4gb
# Maximum run time of job
#SBATCH --time=0-01:00:00
# Configure array parameters
#SBATCH --array 0-999
# Give job a reasonable name
#SBATCH --job-name=fourierT_N30_1e3

# File name for standard output
# (%A will be replaced by the value of SLURM_ARRAY_JOB_ID
# and %a will be replaced by the value of SLURM_ARRAY_TASK_ID)
#SBATCH --output=fourierT_N30_1e3-%A_%a.out
# File name for error output
#SBATCH --error=fourierT_N30_1e3-%A_%a.err

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}


module load phys/qutip/4.6.1

srun python DECODER_q_jumps_iter_w.py ${SLURM_ARRAY_TASK_ID} 0
