#!/bin/bash

# Define the list of arguments
some_ARGS="1 2"
other_ARGS="3 4"
# Iterate over the temperature arguments
for some in $some_ARGS; do
    for other in $other_ARGS; do
        echo "current arguments some=$some, other=$other"   
                ./main.sh $some $other
            done
        done
done
echo "Reached end!"

# Run cleanup code
#echo "Running cleanup code..."
#julia clean_and_store.jl "$(pwd)" 1000
#echo "Cleanup code completed."