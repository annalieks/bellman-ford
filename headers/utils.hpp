#ifndef BELLMAN_FORD_UTILS_HPP
#define BELLMAN_FORD_UTILS_HPP

#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> distr(-10, 1000000);

void generateMatrix(int N, int *&matrix) {
    matrix = (int *) malloc(N * N * sizeof(int));
    for (int i = 0; i < N * N; i++) {
        matrix[i] = distr(gen);
    }
}

#endif //BELLMAN_FORD_UTILS_HPP
