#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

inline double FrobeniusNorm(const std::vector<std::vector<double>> &matrix) {
  double sum = 0.0;
  for (const auto &row : matrix) {
    for (double val : row) {
      sum += val * val;
    }
  }
  return std::sqrt(sum);
}

// returns mean and covariance matrix for selected genes
inline std::pair<std::vector<double>, std::vector<std::vector<double>>>
MeanAndCovariance(const std::vector<double> &selected_genes,
                  size_t num_individuals, size_t d) {
  if (selected_genes.empty() || num_individuals == 0)
    return {};

  std::vector<double> means(d, 0.0);
  for (size_t i = 0; i < num_individuals; ++i) {
    for (size_t j = 0; j < d; ++j) {
      means[j] += selected_genes[i * d + j];
    }
  }

  for (size_t j = 0; j < d; ++j) {
    means[j] /= num_individuals;
  }

  // covariance matrix for fittest individuals
  std::vector<std::vector<double>> cov(d, std::vector<double>(d, 0.0));

  // #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = i; j < d; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < num_individuals; ++k) {
        sum += (selected_genes[k * d + i] - means[i]) *
               (selected_genes[k * d + j] - means[j]);
      }
      cov[i][j] = sum / num_individuals;
      cov[j][i] = cov[i][j];
    }
  }
  return std::make_pair(means, cov);
}

struct RunStats {
  double best;    // BEST  (lowest, since we minimize)
  double worst;   // WORST (highest)
  double mean;    // MEAN
  double median;  // MEDIAN
  double std_dev; // STD DEV
};

inline RunStats ComputeStats(std::vector<double> vals) {
  RunStats s{};
  size_t n = vals.size();

  std::sort(vals.begin(), vals.end());
  s.best = vals.front();
  s.worst = vals.back();

  double sum = 0.0;
  for (double v : vals)
    sum += v;
  s.mean = sum / n;

  s.median = (n % 2 == 0) ? (vals[n / 2 - 1] + vals[n / 2]) / 2.0
                          : vals[n / 2];

  double sum_sq = 0.0;
  for (double v : vals)
    sum_sq += (v - s.mean) * (v - s.mean);
  s.std_dev = std::sqrt(sum_sq / n);

  return s;
}
