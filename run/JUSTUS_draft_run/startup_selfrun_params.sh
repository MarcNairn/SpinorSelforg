#!/bin/bash

# Define the list of arguments
S_ARGS="28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47"
TEMP_ARGS="1 5" #10 12 14 16 18 20"


# Flag to check if there are pending jobs
pending_jobs=false

# Iterate over the temperature arguments
for temp in $TEMP_ARGS; do

    # Count the number of lines in the output (excluding the header)
    num_jobs=$(echo "$(squeue)" | wc -l)
    ((num_jobs--))

    # Print the total number of pending jobs
    echo "Total number of pending jobs: $num_jobs"

    # Check if there are no pending jobs
    if [ "$num_jobs" -eq 0 ]; then
        echo "No pending jobs. Submitting a new job with arguments: temp=$temp..."

        # Run sbatch to submit a new job with the current arguments
        for S in $S_ARGS; do
            sbatch -a 0-49 run/JUSTUS_draft_run/run_justus_static.sh $S $temp 
            echo "New job submitted with arguments: S=$S, temp=$temp."
        done
    else
        echo "There are pending jobs. No action taken for arguments: temp=$temp."
        pending_jobs=true
    fi
done


# Wait for all jobs to finish
echo "Waiting for all jobs to finish..."
while [ $(echo "$(squeue)" | wc -l) -gt 1 ]; do
    sleep 5m  # Adjust the sleep interval as needed
done

echo "Reached end!"

# Run cleanup code
echo "Running cleanup code..."
julia clean_and_store.jl "$(pwd)" 1000
echo "Cleanup code completed."