#include "../include/Generator.h"
#include "../include/Types.h"

#include <random>
#include <vector>

Population RandomPopulation(const Objective &obj, size_t n,
                                   std::mt19937 &rng) {
  auto [lo, hi] = obj.bounds();
  size_t d = obj.dimension();
  std::uniform_real_distribution<double> dis(lo, hi);

  Population population;

  population.all_genes.reserve(n * d);
  population.all_fitnesses.reserve(n);

  population.num_individuals = n;
  population.gene_size = d;

  std::vector<double> cur_gene(d);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      cur_gene[j] = dis(rng);
      population.all_genes.push_back(cur_gene[j]);
    }
    population.all_fitnesses.push_back(obj(cur_gene));
  }
  return population;
}

