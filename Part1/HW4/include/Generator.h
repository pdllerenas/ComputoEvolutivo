#pragma once

#include <array>
#include <cstdint>
#include <random>

#include "Types.h"

Population RandomPopulation(const Objective &obj, size_t n, std::mt19937 &rng);
void GenerateChowLiuPopulation(
    Population &population, const Population &elite, const Objective &obj,
    const DirectedTree &tree, const std::vector<double> &means,
    const std::vector<double> &variances,
    const std::vector<std::vector<double>> &correlations,
    size_t target_population, std::mt19937 &rng);
