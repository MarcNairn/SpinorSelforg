#!/bin/bash
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=2G
# Give job a reasonable name
#SBATCH --job-name=static_pump_Nmc
# File name for standard output (%j will be replaced by job id)
#SBATCH --output=serial_job-%j.out
# File name for error output
#SBATCH --error=serial_job-%j.err

julia run_parallel_justus_static.jl 50 
