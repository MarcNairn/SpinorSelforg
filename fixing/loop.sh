#!/bin/bash


#DEBUGGING SCRIPT TO UNDERSTAND JULIAS ARGUMENT PARSING WHEN THIS GOES FROM ONE BASH SCRIPT TO ANOTHER; LIKE WE DO IN THE MAIN PARAMETER LOOP "startup_selfrun_params" in the run directory"

# D"e"fine t"he list of arguments
some_ARGS="1 2"
other_ARGS="3 4"
# Iterate over the temperature arguments
for some in $some_ARGS; do
    for other in $other_ARGS; do
        echo "current arguments some=$some, other=$other"   
                ./fixing/main.sh "$some" "$other"
            done
        done
done
echo "Reached end!"
