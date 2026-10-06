#pragma once

#include <cstdint>
#include <vector>

#include "Types.h"

std::vector<double> SelectElite(const Population &population, size_t k);
