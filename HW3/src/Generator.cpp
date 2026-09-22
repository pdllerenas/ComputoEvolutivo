#include "../include/Generator.h"

std::array<uint8_t, DATABASE_SIZE>
GenerateDatabase(const std::array<DependencyNode, 6> &probabilities,
                 const std::array<uint16_t, 6> &eval_order, std::mt19937 &rng) {
  std::uniform_real_distribution dis(0.0, 1.0);
  std::array<uint8_t, DATABASE_SIZE> db;

  for (size_t i = 0; i < DATABASE_SIZE; ++i) {
    uint8_t word = 0;
    for (size_t b = 0; b < 6; ++b) {
      int16_t curr_idx = eval_order[b];
      int parent_idx = probabilities[curr_idx].parent_idx;

      double prob = 0.0;
      if (parent_idx == -1) {
        prob = probabilities[curr_idx].p_given_0;
      } else {
        uint8_t parent_bit = (word >> parent_idx) & 1;
        prob = (parent_bit == 1) ? probabilities[curr_idx].p_given_1
                                 : probabilities[curr_idx].p_given_0;
      }
      if (dis(rng) < prob) {
        word |= (1 << curr_idx);
      }
    }
    db[i] = word;
  }
  return db;
}
