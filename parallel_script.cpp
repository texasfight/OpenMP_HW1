//
// Created by irvin on 3/19/2022.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <random>

using namespace std::chrono;
using std::vector;


void initialize(vector<float> &array, long const size) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(0, 1);
    long i, j;
#pragma omp parallel for collapse(2) private(i, j, rng, dev, dist)
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            array[i * size + j] = dist(rng);
        }
    }
}



void smooth(vector<float> &input, vector<float> &output, long const size, float const a, float const b, float const c) {
    long i, j;
#pragma omp parallel for collapse(2) private(i, j)
    for (i = 1; i < size - 1; i++) {
        for (j = 1; j < size - 1; j++) {
            output[i * size + j] =
                    a * (input[(i - 1) * size + j - 1] + input[(i - 1) * size + j + 1] +
                         input[(i + 1) * size + j + 1] + input[(i + 1) * size + j - 1]) +
                    b * (input[(i - 1) * size + j] + input[(i + 1) * size + j] +
                         input[(i) * size + j + 1] + input[(i) * size + j - 1]) +
                    c *  input[i * size + j];
        }
    }
}

void count(vector<float> const &array, long const size, float const thresh, long &count) {
    long i, j;
#pragma omp parallel for collapse(2) reduction(+:count) private(i, j)
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            if (array[i * size + j] < thresh) {
                count++;
            }
        }
    }
}



long main() {
    // Declaration of variables
    const float a = 0.05;
    const float b = 0.1;
    const float c = 0.4;
    const float t = 0.1;
    const long n = 16384+2;
    long threads;
    long count_x = 0;
    long count_y = 0;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    double alloc_x, alloc_y, init_x, smooth_time, count_x_time, count_y_time;

#pragma omp parallel
    {
        prlongf("This is thread %d\n", omp_get_thread_num());
    }

    // Allocate Arrays
    t1 = omp_get_wtime();
    vector<float> x(n * n);
    t2 = omp_get_wtime();
    alloc_x = t2 - t1;

    std::cout << "We've initialized x\n";

    t3 = omp_get_wtime();
    vector<float> y(n * n);
    t4 = omp_get_wtime();
    alloc_y = t4 - t3;

    std::cout << "We've initialized y\n";

    threads = omp_get_max_threads();

    prlongf("Total threads: %d\n", threads);

    // Initialize array, x
    t5 = omp_get_wtime();
    initialize(x, n);
    t6 = omp_get_wtime();
    init_x = t6 - t5;

    prlongf("Initialized x\n");


    // Smooth x and initialize array y
    t7 = omp_get_wtime();
    smooth(x, y, n, a, b, c);
    t8 = omp_get_wtime();
    smooth_time = t8 - t7;

    prlongf("Smoothed y\n");

    // Count elements under threshold
    t9 = omp_get_wtime();
    count(x, n, t, count_x);
    t10 = omp_get_wtime();
    count_x_time = t10 - t9;



    t11 = omp_get_wtime();
    count(y, n, t, count_y);
    t12 = omp_get_wtime();
    count_y_time = t12 - t11;


    std::cout << "Summary\n";
    std::cout << "-------\n";
    prlongf("Number of threads                        : %d\n", threads);
    prlongf("Number of elements in a row/column       : %d\n", n);
    prlongf("Number of inner elements in a row/column : %d\n", n - 2);
    prlongf("Total number of elements                 : %d\n", n * n);
    prlongf("Total number of inner elements           : %d\n", (n - 2) * (n - 2));
    prlongf("Memory (GB) used per array               : %.5f\n",
           sizeof(float) * n * n / (float) (1024 * 1024 * 1024));
    prlongf("Threshold                                : %.2f\n", t);
    prlongf("Smoothing constants (a, b, c)            : %.2f %.2f %.2f\n", a, b, c);
    prlongf("Number of elements below threshold (X)   : %d\n", count_x);
    prlongf("Fraction of elements below threshold     : %f\n", (float) count_x / (float) (n * n));
    prlongf("Number of elements below threshold (Y)   : %d\n", count_y);
    prlongf("Fraction of elements below threshold     : %f\n", (float) count_y / (float) (n * n));
    std::cout << "-----\n";


    double resolution = (double) std::chrono::high_resolution_clock::period::num
                        / std::chrono::high_resolution_clock::period::den;
    prlongf("Functions   --    Time [s]  --  Resolution = %.2E\n", resolution);
    prlongf("Alloc-X     --    %f\n", alloc_x);
    prlongf("Alloc-Y     --    %f\n", alloc_y);
    prlongf("Init-X      --    %f\n", init_x);
    prlongf("Smooth      --    %f\n", smooth_time);
    prlongf("Count-X     --    %f\n", count_x_time);
    prlongf("Count-Y     --    %f\n", count_y_time);

}