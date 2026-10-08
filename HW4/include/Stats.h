#pragma once

#include <vector>

#include "../include/Types.h"

void ComputeMeanVarianceCorrelations(
    const Population &population, std::vector<double> &means,
    std::vector<double> &variances,
    std::vector<std::vector<double>> &correlations);
