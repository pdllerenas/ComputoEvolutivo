#include "../include/Types.h"

#include <vector>

/*
 * Computes the mean and variance of each dimension of population.
 */
void ComputeMeanVarianceCorrelations(const Population &population, std::vector<double> &means,
                       std::vector<double> &variances,
                       std::vector<std::vector<double>> &correlations) {
  size_t d = population.gene_size;
	size_t n = population.num_individuals;

  means.assign(d, 0.0);
 	variances.assign(d, 0.0);
  correlations.assign(d, std::vector<double>(d, 0.0));

  std::vector<double> M2(d, 0.0);
	// co-moments
  std::vector<std::vector<double>> C(d, std::vector<double>(d, 0.0));

	// buffers for Welford update
  std::vector<double> delta(d, 0.0);
  std::vector<double> delta2(d, 0.0);

  // we use this order on the loop to
  // maximize cache efficiency, since
  // each gene is contiguous in memory.
  // going by dimension would certainly
  // generate cache misses.
  for (size_t i = 0; i < n; ++i) {
    double N = static_cast<double>(i + 1);
    for (size_t j = 0; j < d; ++j) {
      double x = population.all_genes[i * d + j];

      delta[j] = x - means[j];
      means[j] += delta[j] / N;
      delta2[j] = x - means[j];
    }
    for (size_t j = 0; j < d; ++j) {
      M2[j] += delta[j] * delta2[j];

      for (size_t k = j + 1; k < d; ++k) {
        C[j][k] += delta[j] * delta2[k];
      }
    }
  }

  if (n > 1) {
    for (size_t j = 0; j < d; ++j) {
      variances[j] = M2[j] / (n - 1);
      correlations[j][j] = 1.0;

      for (size_t k = j + 1; k < d; ++k) {
        if (M2[j] > 1e-12 && M2[k] > 1e-12) {
          double r = C[j][k] / std::sqrt(M2[j] * M2[k]);

          correlations[j][k] = std::clamp(r, -1.0, 1.0);
          correlations[k][j] = correlations[j][k];
        } else {
          correlations[j][k] = 0.0;
          correlations[k][j] = 0.0;
        }
      }
    }
  }
}
