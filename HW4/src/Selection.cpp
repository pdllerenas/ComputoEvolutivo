#include <algorithm>
#include <numeric>
#include <vector>

#include "../include/Selection.h"
#include "../include/Types.h"

Population SelectElite(Population &population, size_t k) {
  size_t n = population.num_individuals;
  size_t d = population.gene_size;

  std::vector<size_t> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  std::partial_sort(indices.begin(), indices.begin() + k, indices.end(),
                    [&population](size_t a, size_t b) {
                      return population.all_fitnesses[a] <
                             population.all_fitnesses[b];
                    });
  std::vector<double> elite_genes, elite_fitnesses;
  elite_genes.reserve(k * d);
  elite_fitnesses.reserve(k);

  for (size_t i = 0; i < k; ++i) {
    size_t og_idx = indices[i];
    auto start_ptr = population.all_genes.begin() + (og_idx * d);
    elite_genes.insert(elite_genes.end(), start_ptr, start_ptr + d);

    elite_fitnesses.push_back(population.all_fitnesses[indices[i]]);
  }
  population.all_genes.clear();
  population.all_fitnesses.clear();

  Population elite;
  elite.all_genes = elite_genes;
  elite.all_fitnesses = elite_fitnesses;
  elite.gene_size = population.gene_size;
  elite.num_individuals = k;
  return elite;
}
