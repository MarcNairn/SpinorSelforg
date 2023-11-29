#!/bin/bash

#SBATCH --job-name=dist-test
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=48
#SBATCH --array 0-3

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export HOME=~

srun julia distributed/parallel_example.jl ${SLURM_CPUS_PER_TASK} ${SLURM_ARRAY_TASK_ID}
