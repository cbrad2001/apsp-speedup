/*
Serial implementation of All-Pairs Shortest Path (Floyd-Warshall) algorithm.
*/

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>



void apspSerial(Graph &g)
{
    // Get number of vertices (n)

    // Create n x n matrices LENGTH and VIA

    // Initialize LENGTH[i][j] for each row i and each col j:
        // If i = j, then LENGTH[i][j] = 0
        // If vertices i and j are neighbors, then LENGTH[i][j] = weight(i,j)
        // Otherwise, LENGTH[i][j] = infinity

    // Initialize VIA[i][j] for each row i and each col j:
        // If i = j, then VIA[i][j] = 0
        // If vertices i and j are neighbors, then VIA[i][j] = weight(i,j)
        // Otherwise, VIA[i][j] = infinity

    // Start serial timer

    // Run apsp serial algorithm:
    // for iteration = 1 to n do
        // for i = 1 to n do
            // for j = 1 to n do
                // if LENGTH[i][iteration] + LENGTH[iteration][j] < LENGTH[i][j] then
                    // LENGTH[i][j] <-- LENGTH[i][iteration] + LENGTH[iteration][j];
                    // VIA[i][j] <-- VIA[i][iteration];

    // Stop timer

    // Output results
}

int main(int argc, char *argv[]) 
{
    cxxopts::Options options(
        "All-Pairs Shortest Path",
        "Calculate all-pairs shortest paths in a graph using serial execution");
    options.add_options(
        "",
        {
            {"nThreads", "Number of Threads",
            cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
            {"inputFile", "Input graph file path",
            cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
        });

    auto cl_options = options.parse(argc, argv);
    uint n_threads = cl_options["nThreads"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

    std::cout << std::fixed;
    std::cout << "Number of Threads : " << n_threads << std::endl;

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    apspSerial(g);

    return 0;
}
