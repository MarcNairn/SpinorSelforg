#!/bin/bash


# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --mem=2G
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
#SBATCH --array 0-0
# Give job a reasonable name
#SBATCH --job-name=static_pump_Nmc
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=serial_job-%j.out
# File name for error output
#SBATCH --error=serial_job-%j.err

srun julia run/JUSTUS_draft_run/run_parallel_justus_static.jl 48 

# you can check the current log of all jobs with the command
# tail -fn 10 [FILENAME]-*
# you can check the job status with
# squeue -i 10

#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export HOME=~
