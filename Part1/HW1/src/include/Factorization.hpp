#pragma once

#include <cmath>
#include <vector>

inline std::vector<std::vector<double>>
Cholesky(const std::vector<std::vector<double>> &A) {
  int N = A.size();
  std::vector<std::vector<double>> L(N, std::vector<double>(N, 0.0));

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < j; ++k) {
        sum += L[i][k] * L[j][k];
      }

      if (i == j) {
        L[i][i] = std::sqrt(std::max(0.0, A[i][i] - sum));
      } else {
        L[i][j] = (1.0 / L[j][j]) * (A[i][j] - sum);
      }
    }
  }
  return L;
}
