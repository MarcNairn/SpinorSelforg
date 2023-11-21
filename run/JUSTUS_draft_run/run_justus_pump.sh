#!/bin/bash

# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=1-00:00:00
#SBATCH --mem=8gb
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
#SBATCH --array 0-32
# Give job a reasonable name
#SBATCH --job-name=pump_range_Nmc
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=pump_range_Nmc-%j.out
# File name for error output
#SBATCH --error=pump_range_Nmc-%j.err

srun julia -p run/JUSTUS_draft_run/run_parallel_justus_pump.jl $SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID
