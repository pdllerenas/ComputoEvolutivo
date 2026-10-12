#pragma once

#include <random>
#include <vector>

#include "Types.h"

namespace ChowLiu {
EDAResult EDA(const Objective &obj, size_t population_size,
              size_t max_generations, std::mt19937 &rng);

} // namespace ChowLiu
