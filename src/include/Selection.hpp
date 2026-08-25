#include <algorithm>
#include <cstdint>
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

std::vector<Individual> Tournament(const std::vector<Individual> &population,
                                   uint16_t match_size = 2,
                                   uint16_t rounds = 1) {
  if (match_size > population.size()) {
    throw std::out_of_range("Match size must not exceed population size.");
  }

  std::vector<Individual> winners(population);
  for (int round = 0; round < rounds; ++round) {
    std::vector<std::vector<Individual>> groups =
        FormGroups(winners, match_size);
    winners = FindWinners(groups);
  }
	return winners;
}
