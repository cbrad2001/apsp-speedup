/*
Parallel implementation of All-Pairs Shortest Path (Floyd-Warshall) algorithm.
*/

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#define DEFAULT_RANDOM_SEED "37"
#define MAX_EDGE_WEIGHT 100

#define IN_TREE 0
#define NOT_IN_TREE 1
#define PIV_LEN 2

static const int INF = 999;

struct processInfo {
	uint start_column;
	uint end_column;
	double process_time_taken;
};

void apspDistributed(Graph &g, uint r_seed, int world_size, int world_rank)
{
    // Initialize timer + time spent
    timer total_timer;
    timer local_timer;
    double total_time_taken = 0.0;
    double local_time_taken = 0.0;
    uintV startNode = 0;
    uintV endNode = 0;

    // Set random generation for given r_seed
    srand(r_seed);

    // Get number of vertices of graph (n)
    uintV n = g.n_;

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

    // std::cout << "Matrices created\n";

    // -------------------------------------------------------------------------------------------
    // Initialize matrices
    for (uintV i = 1; i < n; i++) {
        for (uintV j = 1; j < n; j++) {
            length_curr[i][j] = INF;    // All elements of length_curr and length_next initialized to "infinity"
            length_next[i][j] = INF;
            via_curr[i][j] = INF;     // All elements of via_curr and via_next initialized to 0
            via_next[i][j] = INF;
        }
    }
    for (uintV i = 1; i < n; i++) {
        length_curr[i][i] = 0;      // Same as saying "if i == j then length_curr[i][j] = 0"
        via_curr[i][i] = 0;
        uintE out_degree = g.vertices_[i].getOutDegree();   // Get all outNeighbors (j) of vertex i
        for (uintE deg = 0; deg < out_degree; deg++) {
            uintV j = g.vertices_[i].getOutNeighbor(deg);
            length_curr[i][j] = rand() % MAX_EDGE_WEIGHT + 1;     // Assign "random" edge weight (1 to MAX_EDGE_WEIGHT) to edge(i,j)
            via_curr[i][j] = j;     
        }
    }

    uintV min_nodes = (n - 1) / world_size;
    uintV excess_nodes = (n - 1) % world_size;
    if (world_rank < excess_nodes) {
        startNode = (world_rank * (min_nodes + 1)) + 1;
        endNode = startNode + min_nodes;
    }
    else {
        startNode = (excess_nodes * (min_nodes + 1)) + ((world_rank-excess_nodes) * min_nodes) + 1;
        endNode = startNode + min_nodes - 1;
    }

    printf("Process %d working on nodes %d to %d\n", world_rank, startNode, endNode);

    // local test only (simpleGraph1) - simple 1->4 cyclical graph
    // UNCOMMENT ONLY IF TESTING simpleGraph1
		// length_curr[1][2] = 3;
		// length_curr[2][3] = 5;
		// length_curr[3][4] = 6;
		// length_curr[4][1] = 7;
    
    // for local test only (simpleGraph2) - example graph from https://www.geeksforgeeks.org/floyd-warshall-algorithm-dp-16/
    // UNCOMMENT ONLY IF TESTING simpleGraph2
        // length_curr[1][2] = 4;
        // length_curr[1][4] = 5;
        // length_curr[2][3] = 1;
        // length_curr[2][5] = 6;
        // length_curr[3][1] = 2;
        // length_curr[3][4] = 3;
        // length_curr[4][5] = 2;
        // length_curr[4][3] = 1;
        // length_curr[5][1] = 1;
        // length_curr[5][4] = 4;
    

    // printf("-----------------------------------------\n");
    // printf("initial length[i, j]\n");
    // for (uintV i = 1; i < n; i++) {
    //     for (uintV j = 1; j < n; j++) {
    //         printf("[%3d]", length_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");
    // printf("initial via[i, j]\n");
    // for (uintV i = 1; i < n; i++) {
    //     for (uintV j = 1; j < n; j++) {
    //         printf("[%3d]", via_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");

    // std::cout << "Matrices initialized\n";
    // -------------------------------------------------------------------------------------------
    // Start serial timer
    total_timer.start();
    // -------------------------------------------------------------------------------------------
    // Run apsp serial algorithm:

    // for (uintV iteration = 1; iteration < n; iteration++) {
    //     // Computation phase: Do work on vertices
    //     for (uintV i = 1; i < n; i++) {
    //         for (uintV j = 1; j < n; j++) {
    //             if (length_curr[i][j] > length_curr[i][iteration] + length_curr[iteration][j]
    //             && length_curr[i][iteration] != INF && length_curr[iteration][j] != INF) {
    //                 length_next[i][j] = length_curr[i][iteration] + length_curr[iteration][j];
    //                 via_next[i][j] = via_curr[i][iteration];
    //             }
    //             else {
    //                 length_next[i][j] = length_curr[i][j];
    //                 via_next[i][j] = via_curr[i][j];
    //             }
    //         }
    //     }

    //     // Communication phase: Reset length_next and via_next for next iteration
    //     for (uintV i = 1; i < n; i++) {
    //         for (uintV j = 1; j < n; j++) {
    //             length_curr[i][j] = length_next[i][j];
    //             via_curr[i][j] = via_next[i][j];
    //             length_next[i][j] = INF;
    //             via_next[i][j] = 0;
    //         }
    //     }
    // }

    // Distributed Version
    for (uintV iteration = 1; iteration < n; iteration++) {
        // Computation phase: Do work on vertices
        for (uintV i = startNode; i <= endNode; i++) {
            for (uintV j = 1; j < n; j++) {
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

        // Communication phase: Reset length_next and via_next for next iteration
        for (uintV i = startNode; i <= endNode; i++) {
            for (uintV j = 1; j < n; j++) {
                length_curr[i][j] = length_next[i][j];
                via_curr[i][j] = via_next[i][j];
                length_next[i][j] = INF;
                via_next[i][j] = 0;
            }
        }
    }

    // -------------------------------------------------------------------------------------------
    // Stop timer
    total_time_taken = total_timer.stop();

    // Output results
    if (world_rank == 0){
        std::cout << "thread_id, time_taken" << std::endl;
        std::cout << "0, " << total_time_taken << std::endl;
    }
    
    // printf("-----------------------------------------\n");
    // printf("final length[i, j]\n");
    // for (uintV i = 1; i < n; i++) {
    //     for (uintV j = 1; j < n; j++) {
    //         printf("[%3d]", length_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");
    // printf("final via[i, j]\n");
    // for (uintV i = 1; i < n; i++) {
    //     for (uintV j = 1; j < n; j++) {
    //         printf("[%3d]", via_curr[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("-----------------------------------------\n");
    long long sumLen = 0;
    long long sumVia = 0;
    for (uintV i = 1; i < n; i++) {
        for (uintV j = 1; j < n; j++) {
            sumLen += length_curr[i][j];
        }
    }
    for (uintV i = 1; i < n; i++) {
        for (uintV j = 1; j < n; j++) {
            sumVia += via_curr[i][j];
        }
    }
    printf("Sum Lengths = %lld\n", sumLen);
    printf("Sum Paths = %lld\n", sumVia);

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
            {"inputFile", "Input graph file path",
            cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
            {"rSeed", "Random Seed",
            cxxopts::value<uint>()->default_value(DEFAULT_RANDOM_SEED)}
        });

    auto cl_options = options.parse(argc, argv);
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    uint r_seed = cl_options["rSeed"].as<uint>();

    Graph g;
    // std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    // std::cout << "Created graph\n";

    MPI_Init(NULL, NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0){
        std::cout << "Number of processes : " << world_size << "\n";
        std::cout << std::fixed;
        std::cout << "Random Seed : " << r_seed << "\n";
        std::cout << "rank, start_column, end_column, time_taken\n";
        std::cout << "Number of vertices in graph : " << g.n_ - 1 << "\n";
    }

    if (world_size > g.n_ - 1){
        if (world_rank == 0)
            std::cout << "ERROR: too many processes\n";
        MPI_Finalize();
        return 0;
    }
    apspDistributed(g, r_seed, world_size, world_rank);

    MPI_Finalize();

    return 0;
}

