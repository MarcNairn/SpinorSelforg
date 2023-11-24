#!/bin/bash

# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:01:00
#SBATCH --mem=8gb
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
#SBATCH --array 0-0
# Give job a reasonable name
#SBATCH --job-name=dist-test
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=dist-test-%j.out
# File name for error output
#SBATCH --error=dist-test-%j.err

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export HOME=~

srun julia distributed/parallel_example.jl ${SLURM_CPUS_PER_TASK} ${SLURM_ARRAY_TASK_ID}
