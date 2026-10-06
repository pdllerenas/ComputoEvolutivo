#pragma once

#include <array>
#include <cstdint>
#include <random>

#include "Types.h"

Population RandomPopulation(const Objective &obj, size_t n, std::mt19937 &rng);
