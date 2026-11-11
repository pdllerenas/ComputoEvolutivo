#pragma once

#include <random>
#include <cassert>
#include <vector>

template <typename T> class Population {
private:
  std::vector<T> all_genes;
  std::vector<double> all_fitnesses;
  size_t gene_size = 0;
  size_t num_individuals = 0;

  [[nodiscard]]
  size_t GetFlatIndex(size_t individual_idx) const {
    assert(individual_idx < num_individuals);
    return individual_idx * gene_size;
  }

public:
  Population(size_t num_inds, size_t g_size)
      : num_individuals(num_inds), gene_size(g_size) {
    all_genes.resize(num_individuals * gene_size, T{});
    all_fitnesses.resize(num_individuals, 0.0);
  }

	void IntializeRandom(T lower, T upper, std::mt19937 &rng) {
		std::uniform_int_distribution<T> dis(lower, upper);

		for (size_t i = 0; i < all_genes.size(); ++i) {
			all_genes[i] = dis(rng);
		}
	}

  auto GetGenesBegin(size_t individual_idx) {
    return all_genes.begin() + GetFlatIndex(individual_idx);
  }

  auto GetGenesEnd(size_t individual_idx) {
    return all_genes.begin() + GetFlatIndex(individual_idx) + gene_size;
  }

  auto GetGenesBegin(size_t individual_idx) const {
    return all_genes.cbegin() + GetFlatIndex(individual_idx) + gene_size;
  }

	[[nodiscard]]
  double GetFitness(size_t individual_idx) const {
    assert(individual_idx < num_individuals);
    return all_fitnesses[individual_idx];
  }

  void SetFitness(size_t individual_idx, double fitness) {
    assert(individual_idx < num_individuals);
    all_fitnesses[individual_idx] = fitness;
  }

	[[nodiscard]]
	size_t GetGeneSize() const {
		return gene_size;
	}

	[[nodiscard]]
	size_t GetPopulationSize() const {
		return num_individuals;
	}
};
