#pragma once

#include <vector>

struct Population {
	std::vector<int> all_genes;
	std::vector<double> all_fitnesses;
	size_t gene_size = 0;
	size_t num_individuals = 0;
};
