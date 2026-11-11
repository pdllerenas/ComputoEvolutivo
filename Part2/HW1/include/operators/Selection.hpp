#pragma once

#include <random>

#include "../evolution/Grid1D.hpp"

/*
 * Returns the index to use for crossover operator.
 * It is selected randomly from East-West neighborhood scheme
 */
template <typename T>
size_t BernoulliEastWestSelection(size_t target, size_t west_idx,
                                  size_t east_idx, std::mt19937 &rng) {
  // Randomly Choose Mating Partner
  std::bernoulli_distribution dis(0.5);
  bool random_bool = dis(rng);

  return (random_bool ? west_idx : east_idx);
}

template <typename T>
size_t EastWestSelection(size_t center, size_t west, size_t east,
                         const Grid1D<T> &grid, std::mt19937 &rng,
                         bool is_minimization = true) {

	// Track all neighbors with best fitness for tiebreaker
  std::vector<size_t> best_indices = {center};
  double best_fit = grid.GetReadFitness(center);

  auto evaluate_contender = [&](size_t contender_idx) {
    double contender_fit = grid.GetReadFitness(contender_idx);

		// Assuming we can compare fitnesses directly
    if (contender_fit == best_fit) {
      best_indices.push_back(contender_idx);
    } else {
      bool is_better = is_minimization ? (contender_fit < best_fit)
                                       : (contender_fit > best_fit);
      if (is_better) {
        best_fit = contender_fit;
        best_indices.clear();
        best_indices.push_back(contender_idx);
      }
    }
  };

  evaluate_contender(west);
  evaluate_contender(east);

	// Tiebreaker
  if (best_indices.size() > 1) {
    std::uniform_int_distribution<size_t> dis(0, best_indices.size() - 1);
    return best_indices[dis(rng)];
  }

  return best_indices[0];
}
