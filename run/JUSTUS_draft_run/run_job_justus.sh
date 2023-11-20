#!/bin/bash
for i in {1..31}; do
    echo submitting job $i
		sbatch run_justus_pump.sh $i
		sleep 2s
done
