#pragma once

#include <random>

#include "Factorization.hpp"
#include "Generator.hpp"
#include "Selection.hpp"
#include "Stats.hpp"

inline std::vector<double> EMNA(Population &population, const Objective &obj,
                                std::mt19937 &gen, size_t MAX_EVALS,
                                double sc = 0.3, double TOLERANCE = 1e-4) {
  // while (evals < MAX_EVALS)
  // 	Xs = selection(X, sc)
  // 	[mu, Sigma] = find mean and covariance matrix
  // 	X = generate new population from multivariate normal with parameters mu,
  // Sigma end return best solution

  size_t n = population.num_individuals;
  size_t d = population.gene_size;

  size_t k = static_cast<size_t>(n * sc);

  double optimal = obj.optimal();
  size_t evals = 0;

  bool first_gen = true;
  double initial_norm = 0.0;
  double gamma = 1.0;

  while (evals < MAX_EVALS) {
    std::vector<double> selection = SelectElite(population, k);
    std::vector<double> fittest_individual(selection.begin(),
                                           selection.begin() + d);
    double best_val = obj(fittest_individual);

    std::cerr << evals << "," << best_val << "\n";
    if (std::abs(optimal - best_val) < TOLERANCE) {
      return fittest_individual;
    }
    auto [mu, Sigma] = MeanAndCovariance(selection, k, d);
    double current_norm = FrobeniusNorm(Sigma);
    if (first_gen) {
      initial_norm = current_norm;
      first_gen = false;
    }
    double threshold =
        initial_norm *
        std::pow(static_cast<double>(MAX_EVALS - evals) / MAX_EVALS, gamma);
    if (current_norm > 0.0) {
			// Do theshold * Sigma / ||Sigma|| 
      double scale = threshold / current_norm;
      for (size_t i = 0; i < d; ++i) {
        for (size_t j = 0; j < d; ++j) {
          Sigma[i][j] *= scale;
        }
      }
    }

		// Preserves the condition number of the covariance matrix
    for (size_t i = 0; i < d; ++i) {
    	Sigma[i][i] += 1e-6;
    }

    std::vector<std::vector<double>> L = Cholesky(Sigma);

    MultivariateRandomPopulation(population, obj, fittest_individual, best_val,
                                 L, mu, n, gen);

    evals += n;
  }
  auto it = std::min_element(population.all_fitnesses.begin(),
                             population.all_fitnesses.end());
  size_t best_idx = std::distance(population.all_fitnesses.begin(), it);

  // std::cerr << "EMNA failed to converge in " << MAX_EVALS
  //           << " generations. Best value: " << *it << '\n';
  auto start_ptr = population.all_genes.begin() + (best_idx * d);
  return std::vector<double>(start_ptr, start_ptr + d);
}
