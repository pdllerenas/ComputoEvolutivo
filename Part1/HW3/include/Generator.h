#pragma once

#include <array>
#include <cstdint>
#include <random>

#include "Types.h"

/*
 * @brief Function that generates a database
 * given the dependency lists from empirical observations.
 *
 * We make an optimization on space:
 * we use uint8_t to store the length 6
 * words of the binary alphabet
 *
 * The words generated are of the form
 *
 * 0 0 x6 x5 x4 x3 x2 x1
 */
std::array<uint8_t, DATABASE_SIZE>
GenerateDatabase(const std::array<DependencyNode, 6> &obs,
                 const std::array<uint16_t, 6> &eval_order, std::mt19937 &rng);
