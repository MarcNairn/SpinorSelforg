#!/bin/bash
#SBATCH --job-name=sS-$1
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=0

julia run_parallel_justus_pump.jl 50 $1
