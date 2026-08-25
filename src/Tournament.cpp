#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>

#include "include/Selection.hpp"

int main() {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);

  std::vector<Individual> population;
  population.reserve(100);
  std::cout << "=== Initial population ===\n";
  for (int i = 0; i < 100; ++i) {
    Individual ind = {{0.0 + i, 0.0 + i}, dis(gen)};
    population.push_back(ind);
    std::cout << ind.genes[0] << ", " << ind.fitness << std::endl;
  }
  std::cout << "==========================\n";
	uint8_t match_size = 5;
  std::vector<Individual> winners = Tournament(population, match_size, 2);

  std::cout << "=== Winners ===\n";
  for (auto &winner : winners) {
    std::cout << winner.genes[0] << ", " << winner.fitness << std::endl;
  }
  std::cout << "==========================\n";
  return 0;
}
