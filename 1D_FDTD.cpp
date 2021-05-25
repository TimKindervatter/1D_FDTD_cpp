#include <iostream>
#include "FDTD_engine.h"
#include <chrono>

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    FDTD_engine();
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    double duration_seconds = 1e-9f * duration_ns;
    std::cout << duration_seconds;

    return 0;
}

