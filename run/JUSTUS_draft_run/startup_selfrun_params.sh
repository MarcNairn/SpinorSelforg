#!/bin/bash


# # Define the list of arguments
# S_ARGS="28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47"
# TEMP_ARGS="12 14 16 18" #5 10
# RUN_INDEX="0 1 2 3 4 5 6 7 8 9"

# # Set initial time
# time_allocation="00:30:00"


# for temp in $TEMP_ARGS; do
#     # Iterate over the temperature arguments
#     for run in $RUN_INDEX; do
#         # Iterate over pump strengths
#         for S in $S_ARGS; do
#             # Count the number of lines in the output (excluding the header)
#             num_jobs=$(echo "$(squeue)" | wc -l)
#             ((num_jobs--))
#             # Print the total number of pending jobs
#             echo "Total number of pending jobs: $num_jobs"
#             # Check if there are no pending jobs
#             if [ "$num_jobs" -lt 20 ]; then
#                 echo "Room for new jobs"
#                 # Check the value of S and adjust --time accordingly
#                 if [ "$S" -le 34 ]; then
#                     time_allocation="00:30:00"
#                 elif [ "$S" -le 39 ]; then
#                     # Steeper increase from 30 minutes to 2 and a bit hours
#                     additional_minutes=$((($S - 34) * 20))
#                     time_allocation="00:$(($additional_minutes + 20)):00"
#                 elif [ "$S" -le 44 ]; then
#                     time_allocation="00:240:00"
#                 else
#                     # Finish with a fixed 6 hours for S>44
#                     time_allocation="00:360:00"
#                 fi
#                 sbatch -a 0-49 --time $time_allocation run/JUSTUS_draft_run/run_justus_static.sh $S $temp $run
#                 echo "New job submitted, run=$run, with arguments: S=$S, temp=$temp."
#             else
#                 echo "ERROR: Job queue full, current job will be skipped"
#             fi
#             # Wait for available jobs
#             end_num_jobs = $(echo "$(squeue)" | wc -l)
#             ((num_jobs--))
#             echo "Waiting on queue to free up..."
#             while [end_num_jobs -ge 20 ]; do
#                 sleep 1m  # Adjust the sleep interval as needed
#             done
#         done

#     done
# done

# echo "Reached end!"




## ALTERNATE SCRIPT

# Define the list of arguments
S_ARGS="28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47"
TEMP_ARGS="5 20" # 30 40"
Delta_e=1

# Set initial time
time_allocation="00:30:00"

for temp in $TEMP_ARGS; do
# Iterate over the temperature arguments
    for S in $S_ARGS; do
        num_jobs=$(echo "$(squeue)" | wc -l)
        ((num_jobs--))
        # Print the total number of pending jobs
        echo "Total number of pending jobs: $num_jobs"

        # Check if there are no pending jobs
        if [ "$num_jobs" -eq 0 ]; then
            echo "No pending jobs"

            # Check the value of S and adjust --time accordingly
            if [ "$S" -le 34 ]; then
                time_allocation="00:30:00"
            elif [ "$S" -le 39 ]; then
                additional_minutes=$((($S - 34) * 20))
                time_allocation="00:$(($additional_minutes + 30)):00"
            else
                time_allocation="00:180:00"
            fi

            # Record start time before submitting the job
            start_time=$(date +"%s")

            sbatch -a 0-999 --time $time_allocation run/JUSTUS_draft_run/run_justus_static.sh $S $temp
            echo "New job submitted, with arguments: Delta_e=$Delta_e, S=$S, temp=$temp."

        else
            echo "ERROR: Job queue overran! Current job (S=$S, temp=$temp) will be skipped..."
        fi

	echo "Waiting on jobs to finish..."
        while [ $(echo "$(squeue)" | wc -l) -gt 1 ]; do
            sleep 1s  # Adjust the sleep interval as needed
        done

        # Calculate and print elapsed time for the job
        end_time=$(date +"%s")
        elapsed_time=$((end_time - start_time))
        echo "Elapsed time for S=$S, temp=$temp: $elapsed_time seconds"
        
    done
done

echo "Reached end!"



