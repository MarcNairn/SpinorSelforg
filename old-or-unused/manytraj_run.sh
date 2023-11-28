#!/bin/bash
# Allocate one node
#SBATCH --nodes=1
# Number of program instances to be executed
#SBATCH --ntasks-per-node=1
# 8 GB memory required per node
#SBATCH --mem=4G
# Maximum run time of job
#SBATCH --time=3:00:00
# Give job a reasonable name
#SBATCH --job-name=many_trajectory_draft
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=many_trajectory_draft-%j.out
# File name for error output
#SBATCH --error=many_trajectory_draft-%j.err

# Load software modules as needed, e.g.
module load math/julia/1.7.2

srun julia run/manyrun_draft.jl
# Execute the simulation script in the background
