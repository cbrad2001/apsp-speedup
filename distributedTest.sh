#!/bin/bash

# Check if the input parameter is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <number_of_processes>"
    exit 1
fi

# Assign the input parameter to a variable
num_processes=$1

# Run the command with the specified number of processes
mpirun -n "$num_processes" ./apsp_distributed --inputFile=input_graphs/graph --rSeed=22
