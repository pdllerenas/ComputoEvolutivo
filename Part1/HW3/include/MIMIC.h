#pragma once

#include <array>
#include <vector>

#include "Types.h"

std::array<size_t, WORD_SIZE>
MIMIC(const std::array<uint8_t, DATABASE_SIZE> &db, const EntropyMetrics &em);
