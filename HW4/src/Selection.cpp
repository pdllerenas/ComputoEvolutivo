#include <algorithm>
#include <numeric>
#include <vector>

#include "../include/Selection.h"
#include "../include/Types.h"

std::vector<double> SelectElite(const Population &population, size_t k) {
  size_t n = population.num_individuals;
  size_t d = population.gene_size;

  std::vector<size_t> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  std::partial_sort(indices.begin(), indices.begin() + k, indices.end(),
                    [&population](size_t a, size_t b) {
                      return population.all_fitnesses[a] <
                             population.all_fitnesses[b];
                    });
  std::vector<double> elite_genes;
  elite_genes.reserve(k * d);

  for (size_t i = 0; i < k; ++i) {
    size_t og_idx = indices[i];
    auto start_ptr = population.all_genes.begin() + (og_idx * d);
    elite_genes.insert(elite_genes.end(), start_ptr, start_ptr + d);
  }
  return elite_genes;
}
