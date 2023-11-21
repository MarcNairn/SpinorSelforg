#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --job-name=ex_parallel_simulation
#SBATCH --output=parallel_simulation-%j.out
#SBATCH --error=parallel_simulation-%j.err

srun julia -p $SLURM_CPUS_PER_TASK sampleODE/simplesolve.jl