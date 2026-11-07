#pragma once

#include <chrono>
#include <random>
#include <thread>

class Random {
public:
  Random() {}

  double GetRand(double lower, double upper) {
    // ensure the generator is created once per thread.
    // we do not want race conditions on the in-between step
    // of number generations
    // (https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers)
    static thread_local std::mt19937 generator(
        std::random_device{}() +
        std::hash<std::thread::id>{}(std::this_thread::get_id()));
    std::uniform_real_distribution<double> dis(lower, upper);
    return dis(generator);
  }

  double GetRand() { return GetRand(0.0, 1.0); }
};
