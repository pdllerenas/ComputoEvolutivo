#pragma once

#include <array>
#include <cstdint>
#include <random>

/*
 * @brief Node dependency and probability values
 *
 */
struct DependencyNode {
  int parent_idx;   /**< Index value of parent. -1 for root */
  double p_given_0; /**< Probability of P(X|X_parent = 0). 0 for root */
  double p_given_1; /**< Probability of P(X|X_parent = 1). Marginal for root */
};

/*
 * Dependency nodes for DB1:
 * x1 -> x2 -> x3 -> x4 -> x5 -> x6
 */
constexpr std::array<DependencyNode, 6> CHAIN = {{
    {-1, 0.20, 0.00}, // x1 Root
    {0, 0.92, 0.31},  // x2 Parent is x1
    {1, 0.35, 0.82},  // x3 Parent is x2
    {2, 0.42, 0.15},  // x4 Parent is x3
    {3, 0.28, 0.75},  // x5 Parent is x4
    {4, 0.67, 0.43}   // x6 Parent is x5
}};

// Topological order to evaluate CHAIN
constexpr std::array<uint16_t, 6> CHAIN_EVAL_ORDER = {0, 1, 2, 3, 4, 5};


/*
 * Dependency nodes for DB2:
 * 				x1
 * 			/		\
 * 		 x2		x4
 * 		/ \
 * 	x5  x6
 * 			/
 * 		 x3
 */
constexpr std::array<DependencyNode, 6> TREE = {{
    {-1, 0.20, 0.00}, // x1 Root
    {0, 0.92, 0.31},  // x2 Parent is x1
    {5, 0.35, 0.82},  // x3 Parent is x6
    {0, 0.42, 0.15},  // x4 Parent is x1
    {1, 0.28, 0.75},  // x5 Parent is x2
    {1, 0.67, 0.43}   // x6 Parent is x2
}};

// Topological order to evaluate TREE
constexpr std::array<uint16_t, 6> TREE_EVAL_ORDER = {0, 1, 3, 4, 5, 2};
