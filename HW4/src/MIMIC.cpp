#include <random>

#include "../include/Generator.h"
#include "../include/MIMIC.h"
#include "../include/Selection.h"
#include "../include/Stats.h"
#include "../include/Types.h"

namespace MIMIC {

EDAResult EDA(const Objective &obj, size_t population_size,
              size_t max_generations, std::mt19937 &rng) {
  Population population = RandomPopulation(obj, population_size, rng);

  size_t n = population.num_individuals;
  size_t d = population.gene_size;
  size_t gen = 0;
  bool converged = false;

  std::vector<double> means;
  std::vector<double> variances;
  std::vector<std::vector<double>> correlations;

  Population elite;

  for (; gen < max_generations; ++gen) {
    elite = SelectElite(population, population_size / 2);

    ComputeMeanVarianceCorrelations(elite, means, variances, correlations);

    double max_var = *std::max_element(variances.begin(), variances.end());
    if (max_var < 1e-8) {
      converged = true;
      break;
    }

    std::vector<int> chain;
    chain.reserve(d);
    std::vector<bool> used(d, false);

    // Pick root
    auto iter = std::min_element(variances.begin(), variances.end());
    int root = std::distance(variances.begin(), iter);
    chain.push_back(root);
    used[root] = true;

    while (chain.size() < d) {
      int p = chain.back();
      double best_cond_variance = std::numeric_limits<double>::max();
      int best_idx = -1;
      for (size_t i = 0; i < d; ++i) {
        if (used[i])
          continue;
        double rho = correlations[i][p];
        // H(X|Y) = H(X) - I(X;Y),
        // minimize cond variance: var(x_i | x_p) = var(x_i) * (1-rho^2)
        // is equivalent to minimizing the total entropy
        double cond_var = variances[i] * (1.0 - rho * rho);

        if (cond_var < best_cond_variance) {
          best_cond_variance = cond_var;
          best_idx = i;
        }
      }
      chain.push_back(best_idx);
      used[best_idx] = true;
    }
    DirectedTree dt;
    dt.nodes.resize(d);

    dt.eval_order = chain;

    dt.nodes[chain[0]].parent_idx = -1;

    for (size_t i = 1; i < d; ++i) {
      dt.nodes[chain[i]].parent_idx = chain[i - 1];
    }
    GenerateChowLiuPopulation(population, elite, obj, dt, means, variances,
                              correlations, population_size, rng);
  }
  return {
      std::vector<double>(elite.all_genes.begin(), elite.all_genes.begin() + d),
      elite.all_fitnesses[0], gen, converged};
}
} // namespace MIMIC
