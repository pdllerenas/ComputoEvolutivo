#include <iostream>
#include <random>

#include "../include/ChowLiu.h"
#include "../include/Generator.h"
#include "../include/MIMIC.h"
#include "../include/Stats.h"
#include "../include/Types.h"
#include "../include/Utils.h"

int main() {
  std::random_device rd;
  std::mt19937 rng(rd());

  auto db1 = GenerateDatabase(CHAIN, CHAIN_EVAL_ORDER, rng);
  auto db2 = GenerateDatabase(TREE, TREE_EVAL_ORDER, rng);

  EntropyMetrics em1 = ComputeFrequencies(db1);
  EntropyMetrics em2 = ComputeFrequencies(db2);
  /*
   * Independent Variables
   */
  std::cout << "##### INDEPENDENT #####\n";
  std::cout << "\n=== DATABASE 1 ===\n";

  for (uint8_t i = 0; i < WORD_SIZE; i++) {
    std::cout << "P(X" << i + 1 << " = 1) = "
              << static_cast<double>(em1.marginal_counts[i]) / DATABASE_SIZE
              << '\n';
  }

  double kl_ind_1 = ComputeKLDivergence(db1, em1);
  std::cout << "KL-Divergence = " << kl_ind_1 << "\n";

  std::cout << "\n=== DATABASE 2 ===\n";

  for (uint8_t i = 0; i < WORD_SIZE; i++) {
    std::cout << "P(X" << i + 1 << " = 1) = "
              << static_cast<double>(em2.marginal_counts[i]) / DATABASE_SIZE
              << '\n';
  }

  double kl_ind_2 = ComputeKLDivergence(db2, em2);
  std::cout << "KL-Divergence = " << kl_ind_1 << "\n";

  /*
   * MIMIC
   */
  auto MIMIC_ORDER_1 = MIMIC(db1, em1);

  std::cout << "\n##### MIMIC #####\n";
  std::cout << "\n=== DATABASE 1 ===\n";
  PrintChain(MIMIC_ORDER_1);
  PrintChainProbabilities(MIMIC_ORDER_1, em1);

  double kl_mimic_1 = ComputeKLDivergence(db1, MIMIC_ORDER_1, em1);
  std::cout << "KL-Divergence = " << kl_mimic_1 << "\n";

  std::cout << "\n=== DATABASE 2 ===\n";
  auto MIMIC_ORDER_2 = MIMIC(db2, em2);
  PrintChain(MIMIC_ORDER_2);
  PrintChainProbabilities(MIMIC_ORDER_2, em2);

  double kl_mimic_2 = ComputeKLDivergence(db2, MIMIC_ORDER_2, em2);
  std::cout << "KL-Divergence = " << kl_mimic_2 << "\n";
  /*
   * CHOW-LIU TREE
   */
  std::cout << "\n##### CHOW-LIU #####\n";

  std::cout << "\n=== DATABASE 1 ===\n";
  Graph g1 = BuildChowLiuGraph(db1, em1);
  DirectedTree dt1 = BuildDirectedTree(g1, em1, rng);
  VisualizeTree(dt1);
  double kl_cl_1 = ComputeKLDivergence(db1, dt1);
  std::cout << "KL-Divergence = " << kl_cl_1 << "\n";

  std::cout << "\n=== DATABASE 2 ===\n";
  Graph g2 = BuildChowLiuGraph(db2, em2);
  DirectedTree dt2 = BuildDirectedTree(g2, em2, rng);
  VisualizeTree(dt2);
  double kl_cl_2 = ComputeKLDivergence(db2, dt2);
  std::cout << "KL-Divergence = " << kl_cl_2 << "\n";

  return 0;
}
