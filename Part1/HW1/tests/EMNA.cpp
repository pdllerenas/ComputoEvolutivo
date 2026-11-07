#include <cassert>
#include <iostream>
#include <random>

#include "../src/include/EMNA.hpp"
#include "../src/include/Generator.hpp"
#include "../src/include/Types.hpp"

int main() {
  std::cout << "Testing EMNA\n";

	constexpr size_t NUM_RUNS = 10;
	constexpr int SEED_BASE = 1337;


  size_t d = 2;
  Objective *obj = new Griewank(d);
  size_t n = 1000;

  std::random_device rd;
  std::mt19937 gen(rd());

  const size_t MAX_EVALS = 100'000 * d;

  for (int r = 0; r < NUM_RUNS; ++r) {
    std::mt19937 gen(SEED_BASE + r);
    Population pop = RandomPopulation(*obj, n, gen);
    auto result = EMNA(pop, *obj, gen, MAX_EVALS);
    std::cout << r << "," << ((*obj)(result) - (*obj).optimal()) << "\n";
  }

  // Population population = RandomPopulation(*obj, n, gen);
  // assert(population.gene_size == d);
  // std::vector<double> x = EMNA(population, *obj, gen, MAX_EVALS);

  return 0;
}
