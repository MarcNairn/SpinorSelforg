#!/bin/bash
# Allocate one node
#SBATCH --nodes=1
# Number of program instances to be executed
#SBATCH --ntasks-per-node=1
# 8 GB memory required per node
#SBATCH --mem=4G
# Give job a reasonable name
#SBATCH --job-name=julia_ODE_test
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=serial_job-%j.out
# File name for error output
#SBATCH --error=serial_job-%j.err



srun julia sampleODE_run.jl
# Execute the simulation script in the background
