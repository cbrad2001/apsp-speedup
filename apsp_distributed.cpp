/*
Parallel implementation of All-Pairs Shortest Path (Floyd-Warshall) algorithm.
*/

#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <string.h>

#define DEFAULT_RANDOM_SEED "37"
#define MAX_EDGE_WEIGHT 100

#define IN_TREE 0
#define NOT_IN_TREE 1
#define PIV_LEN 2

static const int INF = 999;

// This helper function returns the rank of the process which has node in its working set
int findDomain(uintV* endNodes, int world_size, int world_rank, uintV node)
{
    for (int i = 0; i < world_size; i++){
        if (node <= endNodes[i]){
            return i;
        }
    }
    return -1;
}

void apspDistributed(Graph &g, uint r_seed, int world_size, int world_rank)
{
    // Initialize timer + time spent
    timer total_timer;
    timer local_timer;
    double total_time_taken = 0.0;
    double local_time_taken = 0.0;

    uintV startNodes[world_size];
    uintV endNodes[world_size];

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
            // TODO: may need to init next vals here?
        }
    }

    uintV min_nodes = n / world_size;
    uintV excess_nodes = n % world_size;
    

    for (int i = 0; i < world_size; i++){
        if (i < excess_nodes) {
            startNodes[i] = i * (min_nodes + 1);
            endNodes[i] = startNodes[i] + min_nodes;
        }
        else {
            startNodes[i] = (excess_nodes * (min_nodes + 1)) + ((i-excess_nodes) * min_nodes);
            endNodes[i] = startNodes[i] + min_nodes - 1;
        }
    }

    // local test only (simpleGraph1) - simple 1->4 cyclical graph
    // UNCOMMENT ONLY IF TESTING simpleGraph1
		// length_curr[0][1] = 3;
		// length_curr[1][2] = 5;
		// length_curr[2][3] = 6;
		// length_curr[3][0] = 7;
    
    // for local test only (simpleGraph2) - example graph from https://www.geeksforgeeks.org/floyd-warshall-algorithm-dp-16/
    // UNCOMMENT ONLY IF TESTING simpleGraph2
        // length_curr[0][1] = 4;
        // length_curr[0][3] = 5;
        // length_curr[1][2] = 1;
        // length_curr[1][4] = 6;
        // length_curr[2][0] = 2;
        // length_curr[2][3] = 3;
        // length_curr[3][4] = 2;
        // length_curr[3][2] = 1;
        // length_curr[4][0] = 1;
        // length_curr[4][3] = 4;

    // -------------------------------------------------------------------------------------------
    // Start timer
    total_timer.start();
    local_timer.start();
    // -------------------------------------------------------------------------------------------
    // Run apsp algorithm:

    for (uintV pivot = 0; pivot < n; pivot++) {    
        
        int pivot_process = findDomain(endNodes, world_size, world_rank, pivot);

        MPI_Bcast(length_curr[pivot], n, MPI_INT32_T, pivot_process, MPI_COMM_WORLD); //propogate the pivot row between all processes

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {
            for (uintV t = 0; t < n; t++){ //perform normal floyd-warshall calculations
                if ((length_curr[i][pivot] + length_curr[pivot][t] < length_curr[i][t])
                && length_curr[i][pivot] != INF && length_curr[pivot][t] != INF){
                    length_next[i][t] = length_curr[i][pivot] + length_curr[pivot][t];
                    via_next[i][t] = via_curr[i][pivot];
                } else {
                    length_next[i][t] = length_curr[i][t];
                    via_next[i][t] = via_curr[i][t];
                }
            }
        }
        
        // Reset length_next and via_next for next pivot
        for (uintV k = startNodes[world_rank]; k <= endNodes[world_rank]; k++) {
            for (uintV j = 0; j < n; j++) {
                length_curr[k][j] = length_next[k][j];
                via_curr[k][j] = via_next[k][j];
                length_next[k][j] = INF;
                via_next[k][j] = INF;
            }
        }

    }

    // -------------------------------------------------------------------------------------------
    // Stop timer
    double timer_stats[world_size] = {0.0};

    local_time_taken = local_timer.stop();


    long long sumLen = 0;
    long long sumVia = 0;

    // UNCOMMENT the following section to view final matrix for len[i, j] and via[i, j]

    // if (world_rank == 0){
    //     for (int i = 0; i < n; i++){
    //         int process_id = findDomain(endNodes, world_size, world_rank, i);
    //         if (process_id != 0){
    //             MPI_Recv(length_curr[i], n, MPI_INT32_T, process_id, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         }
    //     }
    // } else {
    //     for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++) {
    //         MPI_Send(length_curr[i], n, MPI_INT32_T, 0, 44, MPI_COMM_WORLD);
    //     }
    // }

    // if (world_rank == 0){
    //     for (int i = 0; i < n; i++){
    //         int process_id = findDomain(endNodes, world_size, world_rank, i);
    //         if (process_id != 0){
    //             MPI_Recv(via_curr[i], n, MPI_INT32_T, process_id, 34, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         }
    //     }
    // } else {
    //     for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++) {
    //         MPI_Send(via_curr[i], n, MPI_INT32_T, 0, 34, MPI_COMM_WORLD);
    //     }
    // }

    // if (world_rank == 0){
    //     printf("-----------------------------------------\n");
    //     printf("final length[i, j]\n");
    //     for (uintV i = 0; i < n; i++) {
    //         for (uintV j = 0; j < n; j++) {
    //             printf("[%3d]", length_curr[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("-----------------------------------------\n");
    //     printf("final via[i, j]\n");
    //     for (uintV i = 0; i < n; i++) {
    //         for (uintV j = 0; j < n; j++) {
    //             printf("[%3d]", via_curr[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("-----------------------------------------\n");
    // }

    // UNCOMMENT ABOVE

    // Aggregate stats and timers
    for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
        for (uintV j = 0; j < n; j++){
            sumLen += length_curr[i][j];
        }
        for (uintV j = 0; j < n; j++){
            sumVia += via_curr[i][j];
        }
    }
    if (world_rank == 0){
        timer_stats[world_rank] = local_time_taken;
        for (int i = 1; i < world_size; i++){
            long long remoteSumVia;
            long long remoteSumLen;
            MPI_Recv(&remoteSumLen, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&remoteSumVia, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&timer_stats[i], 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sumLen += remoteSumLen;
            sumVia += remoteSumVia;
        }
        
    } else {
        MPI_Send(&sumLen, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&sumVia, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&local_time_taken, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    total_time_taken = total_timer.stop();

    // Print Stats
    if (world_rank == 0){
        std::cout << "rank, start_node, end_node, time_taken\n";
        for (int i = 0; i < world_size; i++){
            std::cout << i << ", " << startNodes[i] << ", " << endNodes[i] << ", " << timer_stats[i] << std::endl;
        }
        
        printf("Sum Lengths = %lld\n", sumLen);
        printf("Sum Paths = %lld\n", sumVia);
        std::cout << "Total Time Taken: " << total_time_taken << std::endl;
    }

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

    MPI_Init(NULL, NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Graph g;
    if (world_rank == 0)
        std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    if (world_rank == 0)
        std::cout << "Created graph\n";

    

    if (world_rank == 0){
        std::cout << "Number of processes : " << world_size << "\n";
        std::cout << std::fixed;
        std::cout << "Random Seed : " << r_seed << "\n";
        std::cout << "Number of vertices in graph : " << g.n_ << "\n";
    }

    if (world_size > g.n_){
        if (world_rank == 0)
            std::cout << "ERROR: too many processes\n";
        MPI_Finalize();
        return 0;
    }
    apspDistributed(g, r_seed, world_size, world_rank);

    MPI_Finalize();

    return 0;
}

