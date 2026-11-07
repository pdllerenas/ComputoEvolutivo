#pragma once

#include <algorithm>
#include <random>
#include <vector>

#include "Types.hpp"

inline Population RandomPopulation(const Objective &obj, size_t n,
                                   std::mt19937 &gen) {
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
    for (size_t j = 0; j < d; ++j) {
      cur_gene[j] = dis(gen);
      population.all_genes.push_back(cur_gene[j]);
    }
    population.all_fitnesses.push_back(obj(cur_gene));
  }

  return population;
}

// Assume L is the result from the Cholesky factorization
// of the covariance matrix
inline void MultivariateRandomPopulation(
    Population &population, const Objective &obj,
    const std::vector<double> &fittest_gene, double fittest_fitness,
    const std::vector<std::vector<double>> &L, const std::vector<double> &means,
    size_t n, std::mt19937 &gen) {
  auto [lo, hi] = obj.bounds();
  size_t d = obj.dimension();
  std::normal_distribution<double> dis(0.0, 1.0);

  population.all_fitnesses.clear();
  population.all_genes.clear();

	population.all_fitnesses.push_back(fittest_fitness);
	population.all_genes.insert(population.all_genes.end(), fittest_gene.begin(), fittest_gene.end());

  std::vector<double> z(d);
  std::vector<double> x(d);

	// only n-1 more individuals, as we kept the fittest
  for (size_t p = 0; p < n - 1; ++p) {
    for (size_t i = 0; i < d; ++i) {
      z[i] = dis(gen);
    }

    for (size_t i = 0; i < d; ++i) {
      x[i] = means[i];
      for (size_t j = 0; j <= i; ++j) {
        x[i] += L[i][j] * z[j];
      }
      x[i] = std::clamp(x[i], lo, hi);
    }
    population.all_fitnesses.push_back(obj(x));
    population.all_genes.insert(population.all_genes.end(), x.begin(), x.end());
  }
}
