#!/bin/bash
# Number of nodes to allocate, always 1
#SBATCH --nodes=1
# Number of MPI instances (ranks) to be executed per node, always 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:10:00
#SBATCH --mem=2G
# Configure array parameters, split job in parts labeled 0-x. (only one job x=0)
#SBATCH --array 0-0
# Give job a reasonable name
#SBATCH --job-name=julia_ODE_test
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=serial_job-%j.out
# File name for error output
#SBATCH --error=serial_job-%j.err



srun julia sampleODE/sampleODE_run.jl
# Execute the simulation script in the background
