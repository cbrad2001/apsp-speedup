/*
Serial implementation of All-Pairs Shortest Path (Floyd-Warshall) algorithm.
*/

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <limits>

static const int INF = std::numeric_limits<int>::max();

void apspSerial(Graph &g)
{
    // Initialize timer + time spent
    timer serial_timer;
    double time_taken = 0.0;

    // Get number of vertices (n)
    uintV n = g.n_;

    // Create n x n matrices LENGTH and VIA
    uintE **length = new uintE[n][n];
    uintE **via = new uintE[n][n];

    // Initialize LENGTH[i][j] and VIA[i][j] for each row i and each col j:
    for (uintV i = 0; i < n; i++) {
        // If i = j, then LENGTH[i][j] = 0, and VIA[i][j] = 0
        // If vertices i and j are neighbors, then LENGTH[i][j] = weight(i,j), and VIA[i][j] = j
        // Otherwise, LENGTH[i][j] = infinity, and VIA[i][j] = infinity
    }

    // Start serial timer
    serial_timer.start();

    // Run apsp serial algorithm:
    // for iteration = 1 to n do
        // for i = 1 to n do
            // for j = 1 to n do
                // if LENGTH[i][iteration] + LENGTH[iteration][j] < LENGTH[i][j] then
                    // LENGTH[i][j] <-- LENGTH[i][iteration] + LENGTH[iteration][j];
                    // VIA[i][j] <-- VIA[i][iteration];

    // Stop timer
    time_taken = serial_timer.stop();

    // Output results
    std::cout << "thread_id, time_taken" << std::endl;
    std::cout << "0, " << time_taken << std::endl;

    // Clean up memory
    delete[][] length;
    delete[][] via;
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
