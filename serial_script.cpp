//
// Created by irvin on 3/19/2022.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <random>

using namespace std::chrono;
using std::vector;


void initialize(vector<float> &array, int const size) {
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



void smooth(vector<float> &input, vector<float> &output, int const size, float const a, float const b, float const c) {
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

void count(vector<float> const &array, int const size, float const thresh, int &count) {
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



int main() {
    // Declaration of variables
    const float a = 0.05;
    const float b = 0.1;
    const float c = 0.4;
    const float t = 0.1;
    const int n = 16384+2;
    int threads;
    int count_x = 0;
    int count_y = 0;
    time_point<system_clock> t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    duration<float> alloc_x, alloc_y, init_x, smooth_time, count_x_time, count_y_time;

    // Allocate Arrays
    t1 = high_resolution_clock::now();
    vector<float> x(n * n);
    t2 = high_resolution_clock::now();
    alloc_x = t2 - t1;


    t3 = high_resolution_clock::now();
    vector<float> y(n * n);
    t4 = high_resolution_clock::now();
    alloc_y = t4 - t3;


    // Initialize array, x
    t5 = high_resolution_clock::now();
    initialize(x, n);
    t6 = high_resolution_clock::now();
    init_x = t6 - t5;


    // Smooth x and initialize array y
    t7 = high_resolution_clock::now();
    smooth(x, y, n, a, b, c);
    t8 = high_resolution_clock::now();
    smooth_time = t8 - t7;


    // Count elements under threshold
    t9 = high_resolution_clock::now();
    count(x, n, t, count_x);
    t10 = high_resolution_clock::now();
    count_x_time = t10 - t9;



    t11 = high_resolution_clock::now();
    count(y, n, t, count_y);
    t12 = high_resolution_clock::now();
    count_y_time = t12 - t11;


    std::cout << "Summary\n";
    std::cout << "-------\n";
    printf("Number of threads                        : %d\n", threads);
    printf("Number of elements in a row/column       : %d\n", n);
    printf("Number of inner elements in a row/column : %d\n", n - 2);
    printf("Total number of elements                 : %d\n", n * n);
    printf("Total number of inner elements           : %d\n", (n - 2) * (n - 2));
    printf("Memory (GB) used per array               : %.5f\n",
           sizeof(float) * n * n / (float) (1024 * 1024 * 1024));
    printf("Threshold                                : %.2f\n", t);
    printf("Smoothing constants (a, b, c)            : %.2f %.2f %.2f\n", a, b, c);
    printf("Number of elements below threshold (X)   : %d\n", count_x);
    printf("Fraction of elements below threshold     : %f\n", (float) count_x / (float) (n * n));
    printf("Number of elements below threshold (Y)   : %d\n", count_y);
    printf("Fraction of elements below threshold     : %f\n", (float) count_y / (float) (n * n));
    std::cout << "-----\n";


    double resolution = (double) std::chrono::high_resolution_clock::period::num
                        / std::chrono::high_resolution_clock::period::den;
    printf("Functions   --    Time [s]  --  Resolution = %.2E\n", resolution);
    printf("Alloc-X     --    %f\n", alloc_x.count());
    printf("Alloc-Y     --    %f\n", alloc_y.count());
    printf("Init-X      --    %f\n", init_x.count());
    printf("Smooth      --    %f\n", smooth_time.count());
    printf("Count-X     --    %f\n", count_x_time.count());
    printf("Count-Y     --    %f\n", count_y_time.count());

}