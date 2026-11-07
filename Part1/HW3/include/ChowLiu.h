#pragma once

#include <random>

#include "Types.h"

Graph BuildChowLiuGraph(const std::array<uint8_t, DATABASE_SIZE> &db,
                 const EntropyMetrics &em);
DirectedTree BuildDirectedTree(Graph &g, const EntropyMetrics &em,
                               std::mt19937 &rng);
