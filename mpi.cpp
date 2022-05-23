#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

#include <mpi.h>
#include "utils.hpp"

using namespace std;

#define INF 1000000

int N; //number of vertices
int p;
int *matrix; // the adjacency matrix
int *dist;
bool hasNegativeCycle = false;

int calculateCoordinate(int x, int y, int n) {
    return x * n + y;
}

int readMatrix(const string &filename) {
    std::ifstream ifStream(filename, std::ifstream::in);
    ifStream >> N;
    matrix = (int *) malloc(N * N * sizeof(int));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            ifStream >> matrix[calculateCoordinate(i, j, N)];
        }
    return 0;
}

void initDistances(int *copyDist, int &copyN) {
    for (int i = 0; i < copyN; i++) {
        copyDist[i] = INF;
    }

    copyDist[0] = 0;
}

int outputResult(const string &file) {
    std::ofstream outStream(file, std::ofstream::out);
    if (!hasNegativeCycle) {
        for (int i = 0; i < N; i++) {
            if (dist[i] > INF)
                dist[i] = INF;
            outStream << dist[i] << '\n';
        }
        outStream.flush();
    } else {
        outStream << "Negative cycle" << endl;
    }
    outStream.close();
    return 0;
}

void relax(int start, int end, int copyN, int *copyMatrix, int *copyDist, bool &rangeChanged) {
    for (int u = start; u < end; u++) {
        for (int v = 0; v < copyN; v++) {
            int weight = copyMatrix[calculateCoordinate(u, v, copyN)];
            if (weight < INF) {
                if (copyDist[u] + weight < copyDist[v]) {
                    copyDist[v] = copyDist[u] + weight;
                    rangeChanged = true;
                }
            }
        }
    }
}

/**
 * Bellman-Ford algorithm. Find the shortest path from vertex 0 to other vertices.
*/
void performBellmanFord(int rank, MPI_Comm comm) {
    int copyN, start, end;
    int *copyMatrix, *copyDist;
    if (rank == 0) {
        copyN = N;
    }

    // broadcast a message to all other processes
    MPI_Bcast(&copyN, 1, MPI_INT, 0, comm);

    // find local task range
    int range = copyN / p;
    start = range * rank;
    end = range * (rank + 1);
    if (rank == p - 1) {
        end = copyN;
    }

    copyMatrix = (int *) malloc(copyN * copyN * sizeof(int));
    copyDist = (int *) malloc(copyN * sizeof(int));

    if (rank == 0) {
        memcpy(copyMatrix, matrix, sizeof(int) * copyN * copyN);
    }
    MPI_Bcast(copyMatrix, copyN * copyN, MPI_INT, 0, comm);

    initDistances(copyDist, N);

    MPI_Barrier(comm);

    bool rangeChanged;
    int iterNum = 0;
    for (int iter = 0; iter < copyN - 1; iter++) {
        rangeChanged = false;
        iterNum++;
        relax(start, end, copyN, copyMatrix, copyDist, rangeChanged);

        MPI_Allreduce(MPI_IN_PLACE, &rangeChanged, 1,
                      MPI_CXX_BOOL, MPI_LOR, comm);
        if (!rangeChanged)
            break;
        MPI_Allreduce(MPI_IN_PLACE, copyDist, copyN,
                      MPI_INT, MPI_MIN, comm);
    }

    if (iterNum == copyN - 1) {
        rangeChanged = false;
        for (int u = start; u < end; u++) {
            for (int v = 0; v < copyN; v++) {
                int weight = copyMatrix[calculateCoordinate(u, v, copyN)];
                if (weight < INF) {
                    if (copyDist[u] + weight < copyDist[v]) {
                        copyDist[v] = copyDist[u] + weight;
                        rangeChanged = true;
                        break;
                    }
                }
            }
        }
        MPI_Allreduce(&rangeChanged, &hasNegativeCycle, 1, MPI_CXX_BOOL, MPI_LOR, comm);
    }

    //step 6: retrieve results back
    if (rank == 0)
        memcpy(dist, copyDist, copyN * sizeof(int));

    //step 7: remember to free memory
    free(copyMatrix);
    free(copyDist);

}

int main(int argc, char **argv) {
    const string inputFile = "../input.txt", outputFile = "../output.txt";
    if (argc <= 1) {
        cerr << "Please, specify number of processes";
        exit(1);
    } else p = stoi(argv[1]);

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int rank; // thread number
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    //only rank 0 process do the I/O
    if (rank == 0) {
        if (argc <= 2) readMatrix(inputFile);
        else generateMatrix(stoi(argv[2]), matrix);
        N = stoi(argv[2]);
        dist = (int *) malloc(sizeof(int) * N);
    }

    //time counter
    double start, end;
    MPI_Barrier(comm);
    start = MPI_Wtime();

    performBellmanFord(rank, comm);
    MPI_Barrier(comm);

    end = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Time: " << (end - start) * 1000 << " ms" << endl;
        outputResult(outputFile);
        free(dist);
        free(matrix);
    }
    MPI_Finalize();
    return 0;
}