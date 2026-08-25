#include <iostream>
#include <cassert>
#include <random>

#include "include/Selection.hpp"

int main() {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);

  std::vector<Individual> population;
  population.reserve(100);
  for (int i = 0; i < 100; ++i) {
    population.push_back({{0.0 + i, 0.0 + i}, dis(gen)});
  }
  std::vector<Individual> winners = Tournament(population);
	assert(winners.size() == 50);
  for (auto &winner : winners) {
    std::cout << winner.genes[0] << ", " << winner.fitness << std::endl;
  }
  return 0;
}
