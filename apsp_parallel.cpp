/*
Parallel implementation of All-Pairs Shortest Path (Floyd-Warshall) algorithm.
*/

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

#define DEFAULT_RANDOM_SEED "37"
#define MAX_EDGE_WEIGHT 100

static const int INF = 999;

void apspParallelTask(Graph &g, uintV **length_curr, uintV **via_curr, uintV **length_next, uintV **via_next, CustomBarrier *barrier, 
                        int thread_id, double *thread_times, uintV start, uintV end)
{
    uintV n = g.n_;
    timer thread_timer;
    thread_timer.start();

    for (uintV iteration = 0; iteration < n; iteration++) {
        // Computation phase: Have thread do work on its vertices
        for (uintV i = start; i <= end; i++) {
            for (uintV j = 0; j < n; j++) {
                if (length_curr[i][j] > length_curr[i][iteration] + length_curr[iteration][j]
                && length_curr[i][iteration] != INF && length_curr[iteration][j] != INF) {
                    length_next[i][j] = length_curr[i][iteration] + length_curr[iteration][j];
                    via_next[i][j] = via_curr[i][iteration];
                }
                else {
                    length_next[i][j] = length_curr[i][j];
                    via_next[i][j] = via_curr[i][j];
                }
            }
        }
        // Have thread wait until all threads complete computation phase
        barrier->wait();

        // Communication phase: Reset length_next and via_next for next iteration
        for (uintV i = start; i <= end; i++) {
            for (uintV j = 0; j < n; j++) {
                length_curr[i][j] = length_next[i][j];
                via_curr[i][j] = via_next[i][j];
                length_next[i][j] = INF;
                via_next[i][j] = INF;
            }
        }
        // Have thread wait until all threads complete communication phase
        barrier->wait();
    }

    thread_times[thread_id] = thread_timer.stop();
}

void apspParallel(Graph &g, int n_threads, uint r_seed)
{
    // Initialize timer + time spent
    timer serial_timer;
    double time_taken = 0.0;

    // Initialize thread data
    std::thread threads[n_threads];
    std::vector<uint> startx(n_threads);
    std::vector<uint> endx(n_threads);
    double *thread_times = new double[n_threads];

    // Set random generation for given r_seed
    srand(r_seed);

    // Get number of vertices of graph (n)
    uintV n = g.n_;

    std::cout << "Number of vertices in graph : " << n << "\n";

    // -------------------------------------------------------------------------------------------
    // Create n x n matrices length_curr, via_curr, length_next, and via_next
    uintV **length_curr = new uintV*[n];
    uintV **via_curr = new uintV*[n];
    uintV **length_next = new uintV*[n];
    uintV **via_next = new uintV*[n];

    for (uintV i = 0; i < n; i++) {
        length_curr[i] = new uintV[n];
        via_curr[i] = new uintV[n];
        length_next[i] = new uintV[n];
        via_next[i] = new uintV[n];
    }

    std::cout << "Matrices created\n";
    // -------------------------------------------------------------------------------------------
    // Initialize matrices
    for (uintV i = 0; i < n; i++) {
        for (uintV j = 0; j < n; j++) {
            length_curr[i][j] = INF;    // All elements of length_curr and length_next initialized to "infinity"
            length_next[i][j] = INF;
            via_curr[i][j] = INF;     // All elements of via_curr and via_next initialized to 0
            via_next[i][j] = INF;
        }
    }
    for (uintV i = 0; i < n; i++) {
        length_curr[i][i] = 0;      // Same as saying "if i == j then length_curr[i][j] = 0"
        via_curr[i][i] = INF;
        uintE out_degree = g.vertices_[i].getOutDegree();   // Get all outNeighbors (j) of vertex i
        for (uintE deg = 0; deg < out_degree; deg++) {
            uintV j = g.vertices_[i].getOutNeighbor(deg);
            length_curr[i][j] = rand() % MAX_EDGE_WEIGHT + 1;     // Assign "random" edge weight (1 to MAX_EDGE_WEIGHT) to edge(i,j)
            via_curr[i][j] = j;     
        }
    }

    std::cout << "Matrices initialized\n";
    // -------------------------------------------------------------------------------------------
    // Initialize barrier
    CustomBarrier *barrier = new CustomBarrier(n_threads);

    // Distribute vertices for each thread task
    uintV min_vertices_per_thread = n / n_threads;
    uintV remaining_vertices = n % n_threads;
    uintV curr_vertex = 0;

    // Initialize thread data arrays (time and number of vertices)
    for (uintV i = 0; i < n_threads; i++) {
        thread_times[i] = 0.0;

        startx[i] = curr_vertex;
        if (remaining_vertices > 0) {
            endx[i] = curr_vertex + min_vertices_per_thread;
            remaining_vertices--;
        }
        else {
            endx[i] = curr_vertex + min_vertices_per_thread - 1;
        }
        curr_vertex = endx[i] + 1;
    }

    // Start serial timer
    serial_timer.start();
    // -------------------------------------------------------------------------------------------
    // Create threads that run apsp on their set of vertices
    for (uintV i = 0; i < n_threads; i++) {
        threads[i] = std::thread(apspParallelTask, std::ref(g), length_curr, via_curr, length_next, via_next, barrier, i, thread_times, startx[i], endx[i]);
    }

    // Wait for threads to join
    for (uintV i = 0; i < n_threads; i++) {
        threads[i].join();
    }
    // -------------------------------------------------------------------------------------------
    // Stop serial timer
    time_taken = serial_timer.stop();

    // UNCOMMENT the following section to view final matrix for len[i, j] and via[i, j]

    // printf("-----------------------------------------\n");
    // printf("final length[i, j]\n");
    // for (uintV i = 0; i < n; i++) {
    //     for (uintV j = 0; j < n; j++) {
    //         printf("[%3d]", length_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");
    // printf("final via[i, j]\n");
    // for (uintV i = 0; i < n; i++) {
    //     for (uintV j = 0; j < n; j++) {
    //         printf("[%3d]", via_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");

    // UNCOMMENT ABOVE

    long long sumLen = 0;
    long long sumVia = 0;
    for (uintV i = 0; i < n; i++) {
        for (uintV j = 0; j < n; j++) {
            sumLen += length_curr[i][j];
        }
    }
    for (uintV i = 0; i < n; i++) {
        for (uintV j = 0; j < n; j++) {
            sumVia += via_curr[i][j];
        }
    }
    printf("Sum Lengths = %lld\n", sumLen);
    printf("Sum Paths = %lld\n", sumVia);

    // Output time for each thread
    std::cout << "thread_id, time_taken" << std::endl;
    for (uintV i = 0; i < n_threads; i++) {
        std::cout << i << ", " << thread_times[i] << "\n";
    }

    // Output total time taken
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    // Clean up memory
    for (uintV i = 0; i < n; i++) {
        delete length_curr[i];
        delete via_curr[i];
        delete length_next[i];
        delete via_next[i];
    }
    delete length_curr;
    delete via_curr;
    delete length_next;
    delete via_next;
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
               "/scratch/input_graphs/graph")},
            {"rSeed", "Random Seed",
            cxxopts::value<uint>()->default_value(DEFAULT_RANDOM_SEED)}
        });

    auto cl_options = options.parse(argc, argv);
    uint n_threads = cl_options["nThreads"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    uint r_seed = cl_options["rSeed"].as<uint>();

    // If command-line parameters are non-positive, print error msg and exit
    if (n_threads <= 0) {
        std::cout << "Error: entered non-positive number of threads and/or max iterations. Shutting down.\n";
        exit(0);
    }

    std::cout << std::fixed;
    std::cout << "Number of Threads : " << n_threads << std::endl;
    std::cout << "Random Seed : " << r_seed << "\n";

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    apspParallel(g, n_threads, r_seed);

    return 0;
}
