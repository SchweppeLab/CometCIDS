#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <ctime>
#include <iostream>
#include <string>

#define BENCH(name, reps) \
    BenchTimer timer(name); \
    for (int64_t i = 0; i < reps; ++i) \

class BenchTimer {
public:
    BenchTimer(std::string name) {
        this->name = name;
        start = clock();
    }

    ~BenchTimer() {
        float seconds = ((float)(clock() - start)) / CLOCKS_PER_SEC;
        std::cout << "benchmark:" << name << ":" << seconds << std::endl;
    }
private:
    std::string name;
    clock_t start;
};

#endif // BENCHMARK_H
