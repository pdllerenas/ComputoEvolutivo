#pragma once

#include <random>
#include <vector>

#include "Types.h"

namespace ChowLiu {
Graph BuildChowLiuGraph(const std::vector<std::vector<double>> &correlations);
DirectedTree BuildDirectedTree(Graph &g, std::mt19937 &rng);
EDAResult EDA(const Objective &obj, size_t population_size,
              size_t max_generations, std::mt19937 &rng);

} // namespace ChowLiu
