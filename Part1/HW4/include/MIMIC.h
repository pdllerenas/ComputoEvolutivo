#pragma once
#include <random>

#include "Types.h"

namespace MIMIC {
EDAResult EDA(const Objective &obj, size_t population_size,
              size_t max_generations, std::mt19937 &rng);
}
