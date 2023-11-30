#!/bin/bash

# Your srun command
JULIA_COMMAND="julia run/JUSTUS_draft_run/run_parallel_justus_static.jl"

# Check if squeue output is empty
if [ -z "$(squeue)" ]; then
    echo "squeue is empty. Running script."

    # Your SLURM script
    sbatch <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=4gb
#SBATCH --array 0-99
#SBATCH --job-name=static_pump_Nmc
#SBATCH --output=static_pump_Nmc-%a_%A.out
#SBATCH --error=static_pump_Nmc-%a_%A.err

$JULIA_COMMAND \$SLURM_ARRAY_TASK_ID 41
EOF

else
    echo "squeue is not empty. Skipping script."
fi