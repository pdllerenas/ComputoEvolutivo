#pragma once

#include "Types.h"
#include <array>
#include <cstdint>
#include <span>

EntropyMetrics ComputeFrequencies(const std::array<uint8_t, DATABASE_SIZE> &db);
double MarginalEntropy(int counts);
double MarginalEntropy(double p);
double ConditionalEntropy(const int *pair_counts);
double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const std::array<size_t, WORD_SIZE> &mimic,
                           const EntropyMetrics &em);
double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const DirectedTree &tree);
double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const EntropyMetrics &em);
