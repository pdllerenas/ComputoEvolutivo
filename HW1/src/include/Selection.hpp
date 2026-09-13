#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "Types.hpp"
// #include "Random.hpp"

inline std::vector<std::vector<Individual>>
FormGroups(std::vector<Individual> &population, uint16_t match_size) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(population.begin(), population.end(), gen);

  std::vector<std::vector<Individual>> groups;
  for (size_t i = 0; i < population.size(); i += match_size) {
    size_t end = std::min(i + match_size, population.size());
    groups.emplace_back(population.begin() + i, population.begin() + end);
  }
  return groups;
}

inline std::vector<Individual>
FindWinners(const std::vector<std::vector<Individual>> &groups) {
  size_t n = groups.size();
  std::vector<Individual> winners;
  winners.reserve(n);
  for (const auto &group : groups) {
    winners.push_back(*std::max_element(group.begin(), group.end()));
  }
  return winners;
}

inline std::vector<Individual>
Tournament(const std::vector<Individual> &population, uint16_t match_size = 2,
           uint16_t rounds = 1) {
  if (match_size > population.size()) {
    throw std::out_of_range("Match size must not exceed population size.");
  }

  std::vector<Individual> winners(population);
  for (uint16_t round = 0; round < rounds; ++round) {
    std::cout << "Round " << round + 1 << std::endl;
    std::vector<std::vector<Individual>> groups =
        FormGroups(winners, match_size);
    winners = FindWinners(groups);
  }
  return winners;
}

// Return a contiguous vector of the k-fittest individuals
inline std::vector<double> SelectElite(const Population &population, size_t k) {
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
