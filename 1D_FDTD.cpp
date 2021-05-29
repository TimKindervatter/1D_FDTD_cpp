#include <iostream>
#include <chrono>
#include <thread>

#include "Eigen\Dense"

#include "FDTD_engine.h"
#include "Utilities.h"

int main()
{
    Eigen::setNbThreads(std::thread::hardware_concurrency());
    std::cout << "Num threads: " << Eigen::nbThreads() << '\n';

    constexpr int SIZE = 10;

    Eigen::Array<double, 1, SIZE> times;

    for (int i = 0; i < SIZE; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();
        FDTD_engine();
        auto stop = std::chrono::high_resolution_clock::now();

        auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
        floating_point_t duration_seconds = 1e-9f * duration_ns;

        times[i] = duration_seconds;

        std::cout << "Run " << i << ": " << duration_seconds << " seconds\n";
    }
        
    std::cout << "Average over " << SIZE << " runs: " << times.mean() << " seconds\n";

    return 0;
}

