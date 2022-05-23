/*
 * This is a openmp version of bellman_ford algorithm
 * Compile: g++ -std=c++11 -o openmp_bellman_ford openmp_bellman_ford.cpp
 * Run: ./openmp_bellman_ford <input file> <number of threads>, you will find the output file 'output.txt'
 * */

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>

#include <omp.h>
#include "utils.hpp"

using namespace std;
using namespace chrono;

#define INF 1000000

int N; //number of vertices
int p = 10; // number of processes
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

void initScope(int *start, int *end, int range) {
#pragma omp parallel for
    for (int i = 0; i < p; i++) {
        start[i] = range * i;
        end[i] = range * (i + 1);
        if (i == p - 1) {
            end[i] = N;
        }
    }
}

void initDistances() {
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        dist[i] = INF;
    }
    dist[0] = 0;
}

void relaxDistances(int *start, int *end, bool *rangeChanged) {
    int threadNum = omp_get_thread_num();
    for (int u = 0; u < N; u++) {
        for (int v = start[threadNum]; v < end[threadNum]; v++) {
            int weight = matrix[calculateCoordinate(u, v, N)];
            if (weight < INF) {
                int newDistance = dist[u] + weight;
                if (newDistance < dist[v]) {
                    rangeChanged[threadNum] = true;
                    dist[v] = newDistance;
                }
            }
        }
    }
}

/**
 * Bellman-Ford algorithm. Find the shortest path from vertex 0 to other vertices.
*/
void performBellmanFord() {
    int start[p], end[p];
    int range = N / p, iterNum = 0;

    omp_set_num_threads(p);
    initScope(start, end, range);
    initDistances();

    bool distanceChanged;
    bool rangeChanged[p];

#pragma omp parallel
    {
        int threadNum = omp_get_thread_num();
        for (int iteration = 0; iteration < N - 1; iteration++) {
            rangeChanged[threadNum] = false;
            relaxDistances(start, end, rangeChanged);
#pragma omp barrier
#pragma omp single
            {
                iterNum++;
                distanceChanged = false;
                for (int t = 0; t < p; t++) {
                    distanceChanged |= rangeChanged[t];
                }
            }
            if (!distanceChanged) {
                break;
            }
        }
    }

    // check negative cycles
    if (iterNum == N - 1) {
        distanceChanged = false;
        for (int u = 0; u < N; u++) {

#pragma omp parallel for reduction(|:distanceChanged)
            for (int v = 0; v < N; v++) {
                int weight = matrix[u * N + v];
                if (weight < INF) {
                    if (dist[u] + weight < dist[v]) {
                        // if we can relax one more step, then we find a negative cycle
                        distanceChanged = true;
                    }
                }
            }
        }
        hasNegativeCycle = distanceChanged;
    }
}

/**
 *
 * @param argc
 * @param argv 1 - number of processes, 2 - number of vertices
 * (if specified, generates random graph)
 * @return
 */
int main(int argc, char **argv) {
    const string inputFile = "../input.txt", outputFile = "../output.txt";
    if (argc <= 1) {
        cerr << "Please, specify number of processes";
        exit(1);
    } else p = stoi(argv[1]);

    if (argc <= 2) readMatrix(inputFile);
    else generateMatrix(stoi(argv[2]), matrix);

    N = stoi(argv[2]);

    dist = (int *) malloc(sizeof(int) * N);

    auto start = steady_clock::now();

    //bellman-ford algorithm
    performBellmanFord();

    auto end = steady_clock::now();
    long duration = duration_cast<nanoseconds>(end - start).count();

    cout << "Time: " << (double) duration / 1000000 << " ms" << endl;
    outputResult(outputFile);
    free(dist);
    free(matrix);
    return 0;
}