#include "../include/Generator.h"
#include "../include/Types.h"

#include <random>
#include <vector>

Population RandomPopulation(const Objective &obj, size_t n, std::mt19937 &rng) {
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
      cur_gene[j] = dis(rng);
      population.all_genes.push_back(cur_gene[j]);
    }
    population.all_fitnesses.push_back(obj(cur_gene));
  }
  return population;
}

void GenerateChowLiuPopulation(
    Population &population, const Population &elite, const Objective &obj,
    const DirectedTree &tree, const std::vector<double> &means,
    const std::vector<double> &variances,
    const std::vector<std::vector<double>> &correlations,
    size_t target_population, std::mt19937 &rng) {

  size_t d = obj.dimension();
  size_t n = target_population;

  population.all_fitnesses.clear();
  population.all_genes.clear();

  population.all_fitnesses.reserve(n);
  population.all_genes.reserve(d * n);

	population.num_individuals = n;
	population.gene_size = d;

  auto [lo, hi] = obj.bounds();

  double fittest_fitness = elite.all_fitnesses[0];
  std::vector<double> fittest_gene(elite.all_genes.begin(),
                                   elite.all_genes.begin() + d);
  // keep the best individual
  population.all_fitnesses.push_back(fittest_fitness);
  population.all_genes.insert(population.all_genes.end(), fittest_gene.begin(),
                              fittest_gene.end());

  std::vector<double> x(d);

  // pre-compute std dev
  std::vector<double> std_devs(d);

  for (size_t i = 0; i < d; ++i) {
    std_devs[i] = std::sqrt(std::max(0.0, variances[i]));
  }

  // create new generation (while keeping best fit)
  for (size_t p = 0; p < n - 1; ++p) {
    // go in topological order
    for (int j : tree.eval_order) {
      int parent = tree.nodes[j].parent_idx;

      // conditional mean and stddev
      double cond_mean = 0.0;
      double cond_stddev = 0.0;

      // for the root, we use the marginal mean and stddev
      if (parent == -1) {
        cond_mean = means[j];
        cond_stddev = std_devs[j];
      } else {
        double rho = correlations[parent][j];
        if (std_devs[parent] > 1e-12) {
          cond_mean = means[j] + rho * (std_devs[j] / std_devs[parent]) *
                                     (x[parent] - means[parent]);
          // 1 - rho^2 gives us the remaining randomness in the children.
          // This shrinks the new stddev, which makes the next iteration
          // move towrads the true mean
          double unexplained_var_ratio = std::max(0.0, 1.0 - (rho * rho));
          cond_stddev = std_devs[j] * std::sqrt(unexplained_var_ratio);
        } else { // when the parent has a near 0 stddev, we also use the
                 // marginals. this occurs when the parent is completely
                 // deterministic, and thus its children do not gain information
                 // from their parent, so they can only rely on their own
                 // estimations
          cond_mean = means[j];
          cond_stddev = std_devs[j];
        }
      }
      // sample from a normal distribution to stay consistent with the
      // maximization of rho^2
      std::normal_distribution<double> dist(cond_mean, cond_stddev);
      x[j] = dist(rng);
      x[j] = std::clamp(x[j], lo, hi);
    }
    population.all_genes.insert(population.all_genes.end(), x.begin(), x.end());
    population.all_fitnesses.push_back(obj(x));
  }
}
