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

struct processInfo {
	uint start_column;
	uint end_column;
	double process_time_taken;
};

struct treeMsg {
    bool inTree;
    uintV pivot;
    uintV targetNode;
};

struct bufferMsg {
    bool isFull;
    uintV* msgBuffer;
};

int findDomain(uintV* endNodes, int world_size, int world_rank, uintV neighbour)
{
    for (int i = 0; i < world_size; i++){
        if (neighbour <= endNodes[i]){
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

    // std::cout << "Matrices created\n";

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

    printf("Process %d working on nodes %d to %d\n", world_rank, startNodes[world_rank], endNodes[world_rank]);

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
    
    // if (world_rank == 0){
    //     printf("-----------------------------------------\n");
    //     printf("initial length[i, j]\n");
    //     for (uintV i = 0; i < n; i++) {
    //         for (uintV j = 0; j < n; j++) {
    //             printf("[%3d]", length_curr[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("-----------------------------------------\n");
    //     printf("initial via[i, j]\n");
    //     for (uintV i = 0; i < n; i++) {
    //         for (uintV j = 0; j < n; j++) {
    //             printf("[%3d]", via_curr[i][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("-----------------------------------------\n");
    // }

    // std::cout << "Matrices initialized\n";
    // -------------------------------------------------------------------------------------------
    // Start serial timer
    total_timer.start();
    // -------------------------------------------------------------------------------------------
    // Run apsp algorithm:

    // Distributed Version (Toueg's Algorithm) - see DC Chapter 5, Page 153
    
    // HEY YOU - PLEASE READ THIS:
    // Prior to testing, run ./input_graphs/SNAPtoBinary input_graphs/simpleGraph1.txt input_graphs/simpleGraph
    // Currently testing for correctness with: mpirun -n 4 ./apsp_distributed --inputFile=input_graphs/simpleGraph --rSeed=20
    // Each process currently only resides over 1 node (marked by i)
    // This is to follow as closely to Toueg's algorithm as possible to ensure correct implementation
    // Once correctly implemented, this should be fairly trivial to implement with node clusters
    // simply will need to use the findDomain() function defined above to find which process to send/rcv to and from
    // CURRENTLY DEADLOCKING AFTER SENDING TREE MESSAGES
    // Avenues to explore:
    // - Firstly, probably best to step through the algorithm with pen and paper but here are some ideas
    // - Are we looping through the correct set of neighbours? ie. maybe need to switch in and out's?
    // - Adding to that, surely we don't need to loop through both in and out neighbours each time if the graph is directed right?
    // - should we update the via matrix when receiving a (in_tree = true) message?
    // - are we correctly manipulating the current and next matrices in the right places?

    for (uintV pivot = 0; pivot < n; pivot++) {         // Line 1
        
        int num_nodes_this_process = endNodes[world_rank] - startNodes[world_rank] + 1;
        std::map<int, MPI_Request*> treeSendsMap;
        std::map<int, MPI_Request*> pivLenSendsMap;
        std::map<int, bool*> isChildMap;
        
        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {
            // Communication Phase: Exchange relevant matrix information between processes

            // Sending whether nodes are part of the sink tree
            uintE in_degree = g.vertices_[i].getInDegree();
            uintE out_degree = g.vertices_[i].getOutDegree();
            MPI_Request* treeSends = new MPI_Request[out_degree];
            for (uintE deg = 0; deg < out_degree; deg++){    // for each inNeighbour
                uintV nbh = g.vertices_[i].getOutNeighbor(deg);  // Line 2
                if (via_curr[i][pivot] == nbh){                 // Line 3
                    treeMsg msg = {
                        true,
                        pivot,
                        nbh
                    };
                    // printf("STEP %d: PROCESS %d SENDING isChild[%d][%d]\n", pivot, world_rank, nbh, i);
                    MPI_Isend(&msg, sizeof(treeMsg), MPI_BYTE, findDomain(endNodes, world_size, world_rank, nbh), nbh, MPI_COMM_WORLD, &treeSends[deg]);     // Line 4
                    //MPI_Send() IN_TREE(PIVOT) to nbh
                } else {
                    treeMsg msg = {
                        false,
                        pivot,
                        nbh
                    };
                    //MPI_Send() NOT_IN_TREE(PIVOT) to nbh
                    MPI_Isend(&msg, sizeof(treeMsg), MPI_BYTE, findDomain(endNodes, world_size, world_rank, nbh), nbh, MPI_COMM_WORLD, &treeSends[deg]);     // Line 5
                }
            }
            treeSendsMap[i] = treeSends;
            
        }

        // printf("Process %d made it past tree sends\n", world_rank);

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
            // bool isChild[n] = {false};
            bool* isChild = new bool[n];
            uintE in_degree = g.vertices_[i].getInDegree();
            
            for (uintE deg = 0; deg < in_degree; deg++){
                uintV nbh = g.vertices_[i].getInNeighbor(deg);
                treeMsg incomingMsg;
                MPI_Recv(&incomingMsg, sizeof(treeMsg), MPI_BYTE, findDomain(endNodes, world_size, world_rank, nbh), i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       // Line 6
                if (incomingMsg.inTree){ //what is the received pivot/nbh used for in the algorithm? is it better to just send a bool?
                        // printf("STEP %d: RECEIVED IN_TREE TRANSMISSION isChild[%d][%d]\n", pivot, i, nbh);
                    isChild[nbh] = true;
                } else {
                    isChild[nbh] = false;
                    // printf("RECEIVED NOT_IN_TREE TRANSMISSION\n");
                }
            }

            isChildMap[i] = isChild;
            
        }

        // printf("Process %d made it past tree recvs\n", world_rank);

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
            //confirm all sends made it
            uintE out_degree = g.vertices_[i].getOutDegree();
            for (uintE deg = 0; deg < out_degree; deg++){
                MPI_Wait(&treeSendsMap[i][deg], MPI_STATUS_IGNORE);       // confirm sends
            }
        }

        // printf("Process %d made it past tree waits\n", world_rank);

        // printf("Process %d has made it past tree sends\n", world_rank);

        // MPI_Request* pivotSends;
        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {

            // MESSAGES WORK FINE UP UNTIL HERE
            // printf("Process %d made it past Tree sends\n", world_rank);

            // Sending relevant rows of the length matrix
            uintE in_degree = g.vertices_[i].getInDegree();
            // MPI_Request* pivLenSends = new MPI_Request[in_degree];
            // MPI_Request pivLenSends[in_degree];
            if (length_curr[i][pivot] != INF){          // Line 7
                if (pivot == i){
                    MPI_Request* pivLenSends = new MPI_Request[in_degree];
                    // pivotSends = new MPI_Request[in_degree];
                    for (uintE deg = 0; deg < in_degree; deg++){       // Line 10
                        uintV nbh = g.vertices_[i].getInNeighbor(deg);
                        if (isChildMap[i][nbh]){     // Line 11
                            int target = findDomain(endNodes, world_size, world_rank, nbh);
                            if (target != world_rank){
                                printf("STEP %d: sending length_curr[%d] P%d -> P%d (node %d) - INIT PIVOT SEND\n", pivot, pivot, world_rank, target, nbh);
                                MPI_Isend(length_curr[pivot], n, MPI_INT32_T, target, 20, MPI_COMM_WORLD, &pivLenSends[deg]);   // Line 12,13,14
                            }
                        }
                    }
                    pivLenSendsMap[i] = pivLenSends;
                }
            }
            

            

        }

        // printf("Process %d made it past piv sends\n", world_rank);

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {
            
            // MESSAGES WORK FINE UP UNTIL HERE
            uintE in_degree = g.vertices_[i].getInDegree();
            if (length_curr[i][pivot] != INF){          // Line 7
                if (pivot != i){                        // Line 8
                    MPI_Request* pivLenSends = new MPI_Request[in_degree];
                    int source = findDomain(endNodes, world_size, world_rank, via_curr[i][pivot]);
                    if (source != world_rank){
                        printf("STEP %d: waiting to receive length_curr[%d] P%d -> P%d (node %d)\n", pivot, pivot, source, world_rank, i);
                        MPI_Recv(length_curr[pivot], n, MPI_INT32_T, source, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        // Line 9
                    }
                    // printf("PROCESS %d RECEIVED LENGTH TRANSMISSION\n", world_rank);
                    for (uintE deg = 0; deg < in_degree; deg++){       // Line 10
                        uintV nbh = g.vertices_[i].getInNeighbor(deg);
                        if (isChildMap[i][nbh]){     // Line 11
                            int target = findDomain(endNodes, world_size, world_rank, nbh);
                            if (target != world_rank){
                                printf("STEP %d: sending length_curr[%d] P%d -> P%d (node %d)\n", pivot, pivot, world_rank, target, nbh);
                                MPI_Isend(length_curr[pivot], n, MPI_INT32_T, target, 20, MPI_COMM_WORLD, &pivLenSends[deg]);   // Line 12,13,14
                            }
                        }
                    }
                    pivLenSendsMap[i] = pivLenSends;
                }
                
            }
            

        }

        // printf("Process %d made it past piv recvs\n", world_rank);

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
            //confirm all sends made it
            uintE in_degree = g.vertices_[i].getInDegree();
            if (length_curr[i][pivot] != INF){
                for (uintE deg = 0; deg < in_degree; deg++){    
                    uintV nbh = g.vertices_[i].getInNeighbor(deg);
                    if (isChildMap[i][nbh] && findDomain(endNodes, world_size, world_rank, nbh) != world_rank){
                        MPI_Wait(&pivLenSendsMap[i][deg], MPI_STATUS_IGNORE); 
                    }
                }
            }
        }


        printf("Process %d made it past piv waits\n", world_rank);

        // MPI_Barrier(MPI_COMM_WORLD);
        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {
            for (uintV t = 0; t < n; t++){          // Lines 15-18
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

        // printf("Process %d updated its arrays\n", world_rank);
        
        // Reset length_next and via_next for next pivot
        for (uintV k = startNodes[world_rank]; k <= endNodes[world_rank]; k++) {
            for (uintV j = 0; j < n; j++) {
                length_curr[k][j] = length_next[k][j];
                via_curr[k][j] = via_next[k][j];
                length_next[k][j] = INF;
                via_next[k][j] = INF;
            }
        }

        // printf("Process %d swapped its arrays\n", world_rank);

        // if (world_rank == 0){
        //     printf("-----------------------------------------\n");
        //     printf("initial length[i, j]\n");
        //     for (uintV i = 0; i < n; i++) {
        //         for (uintV j = 0; j < n; j++) {
        //             printf("[%3d]", length_curr[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("-----------------------------------------\n");
        //     printf("initial via[i, j]\n");
        //     for (uintV i = 0; i < n; i++) {
        //         for (uintV j = 0; j < n; j++) {
        //             printf("[%3d]", via_curr[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("-----------------------------------------\n");
        // }
        // for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
        //     printf("Step %d Row %d Len   : [%3d][%3d][%3d]\n", pivot, i, length_curr[i][0], length_curr[i][1], length_curr[i][2]);
        // }
        // MPI_Barrier(MPI_COMM_WORLD);

        // for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
        //     printf("Step %d Row %d parent: [%3d][%3d][%3d]\n", pivot, i, via_curr[i][0], via_curr[i][1], via_curr[i][2]);
        // }

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++)
        {
            delete treeSendsMap[i];
            delete pivLenSendsMap[i];
            delete isChildMap[i];
        }

        if (world_rank == 0){
            for (int i = 0; i < n; i++){
                int process_id = findDomain(endNodes, world_size, world_rank, i);
                if (process_id != 0){
                    MPI_Recv(length_curr[i], n, MPI_INT32_T, process_id, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        } else {
            for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++) {
                MPI_Send(length_curr[i], n, MPI_INT32_T, 0, 44, MPI_COMM_WORLD);
            }
        }

        if (world_rank == 0){
            for (int i = 0; i < n; i++){
                int process_id = findDomain(endNodes, world_size, world_rank, i);
                if (process_id != 0){
                    MPI_Recv(via_curr[i], n, MPI_INT32_T, process_id, 34, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        } else {
            for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++) {
                MPI_Send(via_curr[i], n, MPI_INT32_T, 0, 34, MPI_COMM_WORLD);
            }
        }

        if (world_rank == 0){
            printf("-----------------------------------------\n");
            printf("Step %d length[i, j]\n", pivot);
            for (uintV i = 0; i < n; i++) {
                for (uintV j = 0; j < n; j++) {
                    printf("[%3d]", length_curr[i][j]);
                }
                printf("(%d)\n", findDomain(endNodes, world_size, world_rank, i));
            }
            printf("-----------------------------------------\n");
            printf("Step %d via[i, j]\n", pivot);
            for (uintV i = 0; i < n; i++) {
                for (uintV j = 0; j < n; j++) {
                    printf("[%3d]", via_curr[i][j]);
                }
                printf("(%d)\n", findDomain(endNodes, world_size, world_rank, i));
            }
            printf("-----------------------------------------\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    // printf("all done!\n");

    // -------------------------------------------------------------------------------------------
    // Stop timer
    total_time_taken = total_timer.stop();

    // Output results
    if (world_rank == 0){
        // std::cout << "thread_id, time_taken" << std::endl;
        // std::cout << "0, " << total_time_taken << std::endl;
    }
    // printf("Process %d final Len: [%3d][%3d][%3d][%3d]\n", world_rank, length_curr[world_rank][0], length_curr[world_rank][1], length_curr[world_rank][2], length_curr[world_rank][3]);
    
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
    // long long sumLen = 0;
    // long long sumVia = 0;
    // for (uintV i = 0; i < n; i++) {
    //     for (uintV j = 0; j < n; j++) {
    //         sumLen += length_curr[i][j];
    //     }
    // }
    // for (uintV i = 0; i < n; i++) {
    //     for (uintV j = 0; j < n; j++) {
    //         sumVia += via_curr[i][j];
    //     }
    // }
    // printf("Sum Lengths = %lld\n", sumLen);
    // printf("Sum Paths = %lld\n", sumVia);
    long long sumLen = 0;
    long long sumVia = 0;

    for (int i = startNodes[world_rank]; i <= endNodes[world_rank]; i++){
        for (uintV j = 0; j < n; j++){
            sumLen += length_curr[i][j];
        }
        for (uintV j = 0; j < n; j++){
            sumVia += via_curr[i][j];
        }
    }
    if (world_rank == 0){
        std::cout << "thread_id, time_taken" << std::endl;
        std::cout << "0, " << total_time_taken << std::endl;
        for (int i = 1; i < world_size; i++){
            long long remoteSumLen;
            MPI_Recv(&remoteSumLen, 1, MPI_LONG_LONG, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sumLen += remoteSumLen;
        }
        for (int i = 1; i < world_size; i++){
            long long remoteSumVia;
            MPI_Recv(&remoteSumVia, 1, MPI_LONG_LONG, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sumVia += remoteSumVia;
        }
        printf("Sum Lengths = %lld\n", sumLen);
        printf("Sum Paths = %lld\n", sumVia);
    } else {
        MPI_Send(&sumLen, 1, MPI_LONG_LONG, 0, 99, MPI_COMM_WORLD);
        MPI_Send(&sumVia, 1, MPI_LONG_LONG, 0, 99, MPI_COMM_WORLD);
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

