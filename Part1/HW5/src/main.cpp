#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "../include/Takahama.h"
#include "../include/Types.h"

void print_ranking(size_t gen, double eps,
                   const std::vector<size_t> &sorted_indices,
                   const Population &pop) {
  std::cout << "\nGeneration " << gen << " | Epsilon Threshold: " << eps
            << "\n";
  std::cout << "Rank | (x, y)          | f(x)       | v(x)\n";
  std::cout << "--------------------------------------------------\n";

  for (size_t rank = 0; rank < sorted_indices.size(); ++rank) {
    size_t idx = sorted_indices[rank];
    double x = pop.all_genes[idx * pop.gene_size];
    double y = pop.all_genes[idx * pop.gene_size + 1];

    std::cout << std::setw(4) << rank + 1 << " | "
              << "(" << std::setw(5) << x << ", " << std::setw(4) << y << ") | "
              << std::setw(10) << pop.all_fitnesses[idx] << " | "
              << std::setw(10) << pop.all_violations[idx] << "\n";
  }
}

int main() {
  Himmelblau himmelblau_obj(2);
  HimmelblauConstraints himmelblau_constraints;

  // 5 non-feasible points within the search space (x >= 0, y >= 0)
  std::vector<std::vector<double>> test_points = {
      {10.0, 5.0}, {6.0, 6.0}, {0.0, 25.0}, {8.0, 1.0}, {5.0, 10.0}};

  // Init population
  Population pop;
  pop.num_individuals = test_points.size();
  pop.gene_size = 2;
  pop.all_genes.reserve(pop.num_individuals * pop.gene_size);

  for (const auto &pt : test_points) {
    pop.all_genes.push_back(pt[0]);
    pop.all_genes.push_back(pt[1]);
  }

  batch_evaluate(pop, himmelblau_obj, himmelblau_constraints);

  double eps0 = get_initial_eps0(pop.all_violations);
  size_t max_generations = 1000;
  size_t T_control = 800;
  double cp = 4.0;

  Takahama schedule(eps0, T_control, cp);

  std::vector<size_t> indices(pop.num_individuals);

  std::vector<size_t> test_generations = {0, T_control / 2, T_control};

  std::cout << std::fixed << std::setprecision(2);

  for (size_t t : test_generations) {
    double current_eps = schedule(t);

    std::iota(indices.begin(), indices.end(), 0);

    std::sort(
        indices.begin(), indices.end(),
        EpsilonCompare(pop.all_fitnesses, pop.all_violations, current_eps));

    print_ranking(t, current_eps, indices, pop);
  }

  return 0;
}
