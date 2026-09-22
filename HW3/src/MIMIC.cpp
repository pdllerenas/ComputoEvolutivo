#include <iostream>

#include "../include/MIMIC.h"
#include "../include/Stats.h"

std::array<size_t, WORD_SIZE>
MIMIC(const std::array<uint8_t, DATABASE_SIZE> &db, const EntropyMetrics &em) {

  std::array<size_t, WORD_SIZE> MIMIC_ORDER;
  std::vector<uint8_t> NODES = {0, 1, 2, 3, 4, 5};

  // Find lowest entropy node, this way we start with the
  // node which has found the best 'patterns'. This defines the root node
  double min_entropy = std::numeric_limits<double>::max();
  uint8_t root_node = 0;

  for (size_t i = 0; i < WORD_SIZE; ++i) {
    double p = static_cast<double>(em.marginal_counts[i]) / DATABASE_SIZE;
    double h = MarginalEntropy(p);
    if (h < min_entropy) {
      min_entropy = h;
      root_node = i;
    }
  }

  MIMIC_ORDER[0] = root_node;
  std::erase(NODES, root_node);

  for (int i = 1; i < WORD_SIZE; i++) {
    uint8_t parent_node = MIMIC_ORDER[i - 1];

    double best_conditional_entropy = std::numeric_limits<double>::max();
    uint8_t best_node = 0;

    for (uint8_t candidate : NODES) {
      // flat array index for conditional entropy
      int offset = ((candidate * WORD_SIZE) + parent_node) * 4;

      double cond_entropy = ConditionalEntropy(&em.joint_counts[offset]);
      if (cond_entropy < best_conditional_entropy) {
        best_conditional_entropy = cond_entropy;
        best_node = candidate;
      }
    }
    MIMIC_ORDER[i] = best_node;
    std::erase(NODES, best_node);
  }

	return MIMIC_ORDER;
}
