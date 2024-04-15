### Group 7 - 431 Project: APSP - Floy-Warshall Algorithm

## Compilation Instruction:

Unzip submission

`tar xvzf proj.tar.gz`

Run make to compile executables

`make`

## Graph Generation

The input_graphs/ directory contains 4 .txt files which can be used for generating graphs

    10Nodes.txt
    100Nodes.txt
    1000Nodes.txt
    5000Nodes.txt

To generate a graph from the txt files, run the following command

`./input_graphs/SNAPtoBinary input_graphs/{$graph_txt_file} input_graphs/{graph_name}`

For example:

`./input_graphs/SNAPtoBinary input_graphs/10Nodes.txt input_graphs/10Nodes`

Feel free to also generate your own graph using our python script

`cd input_graphs/`

`python3 graphgen.py 50 100`

`cd ..`

`./input_graphs/SNAPtoBinary input_graphs/testGraph.txt input_graphs/testGraph`

This will create a cutsome graph with 50 nodes and 100 edges

## Run instructions

1. To run the serial algorithm, run:

`./apsp_serial  --inputFile input_graphs/{graph_name} --rSeed {$seed}`

For Example
`./apsp_serial  --inputFile input_graphs/10Nodes --rSeed 22`

2. To run the parallel threaded algorithm, run:

`./apsp_parallel --inputFile input_graphs/{graph_name} --rSeed {$seed} --nThreads {$n_threads}`

For Example
`./apsp_parallel --inputFile input_graphs/10Nodes --rSeed 22 --nThreads 8`

3. To run the distributed threaded algorithm, run:

`mpirun -n {n_processes} ./apsp_distributed --inputFile=input_graphs/{graph_name} --rSeed={$seed}`

For Example
`mpirun -n 8 ./apsp_distributed --inputFile input_graphs/10Nodes --rSeed 22`

## Important Note about input

Ensure that you use the same seed when comparing outputs between different versions of the algorithms. Edge weights are assigned to the graphs at runtime based on the seed. Therefore to accurately compare output and runtimes between different versions of the algorithm, the same seed must be used to generate the same input.

## Output

Each version of the algorithm will print the sum of the length and via matrices to prove correctness of the algorithm between versions, as well as runtimes for each thread/process and total runtime. Initially, we also printed the matrices to stdout however this did not scale well for large number of processes.

The option to print out matrices to stdout still exists if the user wishes to enable it, they must simply uncomment the sections in each algorithm file which pertains to printing matrix output. The code blocks are found in the files as follows:

Serial - Line 126

Parallel - Line 173

Distributed - Line 171

