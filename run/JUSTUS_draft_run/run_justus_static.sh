#!/bin/bash


# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4gb
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
#SBATCH --array 0-99
# Give job a reasonable name
#SBATCH --job-name=static_pump_Nmc
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=static_pump_Nmc-%a_%A.out
# File name for error output
#SBATCH --error=static_pump_Nmc-%a_%A.err

srun julia run/JUSTUS_draft_run/run_parallel_justus_static.jl ${SLURM_ARRAY_TASK_ID} 40

#ARGS[2] represents the pumping strength

# you can check the current log of all jobs with the command
# tail -fn 10 [FILENAME]-*
# you can check the job status with
# squeue -i 10

#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export HOME=~
