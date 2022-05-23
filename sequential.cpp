#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>

#include "utils.hpp"

using namespace std;
using namespace chrono;

#define INF 1000000

int N; //number of vertices
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

void initDistances() {
    for (int i = 0; i < N; i++) {
        dist[i] = INF;
    }
    dist[0] = 0;
}

void relaxDistances(bool& distanceChanged) {
    for (int u = 0; u < N; u++) {
        for (int v = 0; v < N; v++) {
            int weight = matrix[calculateCoordinate(u, v, N)];
            if (weight < INF) {
                int newDistance = dist[u] + weight;
                if (newDistance < dist[v]) {
                    distanceChanged = true;
                    dist[v] = newDistance;
                }
            }
        }
    }
}

void checkNegativeCycle() {
    for (int u = 0; u < N; u++) {
        for (int v = 0; v < N; v++) {
            int weight = matrix[calculateCoordinate(u, v, N)];
            if (weight < INF) {
                if (dist[u] + weight < dist[v]) {
                    // if we can relax one more step, then we find a negative cycle
                    hasNegativeCycle = true;
                    return;
                }
            }
        }
    }
}

/**
 * Bellman-Ford algorithm. Find the shortest path from vertex 0 to other vertices.
*/
void performBellmanFord() {
    initDistances();

    bool distanceChanged;
    for (int i = 0; i < N - 1; i++) {// n - 1 iteration
        distanceChanged = false;
        relaxDistances(distanceChanged);
        if (!distanceChanged) {
            return;
        }
    }

    //do one more iteration to check negative cycles
    checkNegativeCycle();
}

int main(int argc, char **argv) {
    const string inputFile = "../input.txt", outputFile = "../output.txt";

    if (argc <= 1) readMatrix(inputFile);
    else generateMatrix(stoi(argv[1]), matrix);

    N = stoi(argv[1]);

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