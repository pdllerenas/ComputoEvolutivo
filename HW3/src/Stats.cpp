#include <cmath>

#include "../include/Stats.h"
#include "../include/Types.h"
#include "../include/Utils.h"

static void UpdateMarginals(uint8_t word, std::vector<int> &marginal_counts) {
  for (size_t i = 0; i < WORD_SIZE; ++i) {
    if ((word >> i) & 1) {
      marginal_counts[i]++;
    }
  }
}

static void UpdateJoints(uint8_t word, std::vector<int> &joint_counts) {
  for (size_t i = 0; i < WORD_SIZE; ++i) {
    for (size_t j = 0; j < WORD_SIZE; ++j) {
      uint8_t bit_i = (word >> i) & 1;
      uint8_t bit_j = (word >> j) & 1;
      int state_idx =
          (bit_i << 1) | bit_j; // gives some number in [0,1,2,3]: 00 -> 0, 01
                                // -> 1, 10 -> 2, 11 -> 3
      int flat_idx = ((i * WORD_SIZE) + j) * 4 + state_idx;
      joint_counts[flat_idx]++;
    }
  }
}

EntropyMetrics
ComputeFrequencies(const std::array<uint8_t, DATABASE_SIZE> &db) {
  EntropyMetrics metrics;
  metrics.marginal_counts.assign(WORD_SIZE, 0);
  metrics.joint_counts.assign(4 * WORD_SIZE * WORD_SIZE, 0);

  for (uint8_t word : db) {
    UpdateMarginals(word, metrics.marginal_counts);
    UpdateJoints(word, metrics.joint_counts);
  }
  return metrics;
}

/*
 * @brief Calcultes the marginal entropy given the count of a variable
 * (Assuming binary random variable)
 */
double MarginalEntropy(int counts) {
  double p = static_cast<double>(counts) / DATABASE_SIZE;
  return MarginalEntropy(p);
}

/*
 * @brief Calcultes the marginal entropy given the marginal probability
 * (Assuming binary random variable)
 */
double MarginalEntropy(double p) {
  if (p <= 1e-12 || p + 1e-12 >= 1)
    return 0;
  return -(p * std::log2(p) + (1 - p) * std::log2(1 - p));
}

/*
 * @brief Calcultes the marginal entropy given the marginal probability
 * (Assuming binary random variable)
 */
double ConditionalEntropy(const int *pair_counts) {
  int total_N =
      pair_counts[0] + pair_counts[1] + pair_counts[2] + pair_counts[3];
  if (total_N == 0)
    return 0.0;

  double p00 = static_cast<double>(pair_counts[0]) / total_N;
  double p01 = static_cast<double>(pair_counts[1]) / total_N;
  double p10 = static_cast<double>(pair_counts[2]) / total_N;
  double p11 = static_cast<double>(pair_counts[3]) / total_N;

  double p_j0 = p00 + p10;
  double p_j1 = p01 + p11;

  double entropy = 0.0;

  if (p_j0 > 1e-12) {
    double p_i0_given_j0 = p00 / p_j0;
    double p_i1_given_j0 = p10 / p_j0;

    double h_cond_0 = 0.0;
    if (p_i0_given_j0 > 1e-12)
      h_cond_0 -= p_i0_given_j0 * std::log2(p_i0_given_j0);
    if (p_i1_given_j0 > 1e-12)
      h_cond_0 -= p_i1_given_j0 * std::log2(p_i1_given_j0);

    entropy += p_j0 * h_cond_0;
  }

  if (p_j1 > 1e-12) {
    double p_i0_given_j1 = p01 / p_j1;
    double p_i1_given_j1 = p11 / p_j1;

    double h_cond_1 = 0.0;
    if (p_i0_given_j1 > 1e-12)
      h_cond_1 -= p_i0_given_j1 * std::log2(p_i0_given_j1);
    if (p_i1_given_j1 > 1e-12)
      h_cond_1 -= p_i1_given_j1 * std::log2(p_i1_given_j1);

    entropy += p_j1 * h_cond_1;
  }

  return entropy;
}

double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const DirectedTree &tree) {
  std::array<int, 64> state_counts = {0};
  for (uint8_t word : db) {
    state_counts[word]++;
  }

  double kl_div = 0.0;

  for (int x = 0; x < 64; ++x) {
    if (state_counts[x] == 0)
      continue;

    double p_x = static_cast<double>(state_counts[x]) / DATABASE_SIZE;
    double q_x = 1.0;
    for (int i = 0; i < WORD_SIZE; ++i) {
      int bit_i = (x >> i) & 1;
      int parent = tree.nodes[i].parent_idx;

      if (parent == -1) {
        q_x *= (bit_i == 1) ? tree.nodes[i].p_given_1
                            : (1.0 - tree.nodes[i].p_given_1);
      } else {
        int bit_parent = (x >> parent) & 1;
        double p1 = (bit_parent == 1) ? tree.nodes[i].p_given_1
                                      : tree.nodes[i].p_given_0;
        q_x *= (bit_i == 1) ? p1 : (1.0 - p1);
      }
    }
    if (q_x > 1e-12) {
      kl_div += p_x * std::log2(p_x / q_x);
    }
  }
  return kl_div;
}

double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const std::array<size_t, WORD_SIZE> &mimic,
                           const EntropyMetrics &em) {
  DirectedTree dt = MimicToTree(mimic, em);
  return ComputeKLDivergence(db, dt);
}

double ComputeKLDivergence(const std::array<uint8_t, DATABASE_SIZE> &db,
                           const EntropyMetrics &em) {
  DirectedTree dt = IndependentToTree(em);
  return ComputeKLDivergence(db, dt);
}
