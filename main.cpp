#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>

using namespace std::chrono;
using std::vector;


void initialize(vector<float> &array, int const size) {
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            array[i * size + j] = rand() / (float) RAND_MAX;
        }
    }
}



void smooth(vector<float> &input, vector<float> &output, int const size, float const a, float const b, float const c) {
    for (int i = 1; i < size - 1; i++) {
        for (int j = 1; j < size - 1; j++) {
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
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
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

    // Allocate Arrays
    auto t1 = high_resolution_clock::now();
    vector<float> x(n*n);
    auto t2 = high_resolution_clock::now();;
    duration<float> alloc_x = t2 - t1;


    t1 = high_resolution_clock::now();
    vector<float> y(n*n);
    t2 = high_resolution_clock::now();
    duration<float> alloc_y = t2 - t1;


    // Initialize array, x
    t1 = high_resolution_clock::now();
    initialize(x, n);
    t2 = high_resolution_clock::now();
    duration<float> init_x = t2 - t1;

    // Smooth x and initialize array y
    t1 = high_resolution_clock::now();
    smooth(x, y, n, a, b, c);
    t2 = high_resolution_clock::now();
    duration<float> smooth = t2 - t1;

    // Count elements under threshold
    t1 = high_resolution_clock::now();
    int count_x = 0;
    count(x, n, t, count_x);
    t2 = high_resolution_clock::now();
    duration<float> count_x_time = t2 - t1;


    t1 = high_resolution_clock::now();
    int count_y = 0;
    count(y, n, t, count_y);
    t2 = high_resolution_clock::now();
    duration<float> count_y_time = t2 - t1;

    std::cout << "Summary\n";
    std::cout << "-------\n";
    printf("Number of elements in a row/column       : %d\n", n);
    printf("Number of inner elements in a row/column : %d\n", n - 2);
    printf("Total number of elements                 : %d\n", n * n);
    printf("Total number of inner elements           : %d\n", (n - 2) * (n - 2));
    printf("Memory (GB) used per array               : %.5f\n", sizeof(float) * n * n / (float) (1024 * 1024 * 1024));
    printf("Threshold                                : %.2f\n", t);
    printf("Smoothing constants (a, b, c)            : %.2f %.2f %.2f\n", a, b, c);
    printf("Number of elements below threshold (X)   : %d\n", count_x);
    printf("Fraction of elements below threshold     : %f\n",  (float) count_x / (float) (n * n));
    printf("Number of elements below threshold (Y)   : %d\n", count_y);
    printf("Fraction of elements below threshold     : %f\n", (float) count_y / (float) (n * n));
    std::cout << "-----\n";


    double resolution = (double) std::chrono::high_resolution_clock::period::num
                        / std::chrono::high_resolution_clock::period::den;
    printf("Functions   --    Time [s]  --  Resolution = %.2E\n", resolution);
    printf("Alloc-X     --    %f\n", alloc_x.count());
    printf("Alloc-Y     --    %f\n", alloc_y.count());
    printf("Init-X      --    %f\n", init_x.count());
    printf("Smooth      --    %f\n", smooth.count());
    printf("Count-X     --    %f\n", count_x_time.count());
    printf("Count-Y     --    %f\n", count_y_time.count());

}