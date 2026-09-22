#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <string>

#include "ChowLiu.h"
#include "Types.h"

inline void Print_uint8_Bits(uint8_t n) {
  for (int i = 7; i >= 0; --i) {
    std::cout << ((n >> i) & 1);
  }
  std::cout << '\n';
}

inline void PrintChain(const std::array<size_t, WORD_SIZE> &MIMIC_ORDER) {
  for (int i = 0; i < MIMIC_ORDER.size(); ++i) {
    std::cout << MIMIC_ORDER[i] + 1;
    if (i < MIMIC_ORDER.size() - 1) {
      std::cout << " -> ";
    }
  }
  std::cout << '\n';
}
// Forward declaration for the recursive helper
static void
PrintTreeRecursive(int node, const std::vector<std::vector<int>> &children,
                   const std::array<DependencyNode, WORD_SIZE> &nodes,
                   std::string prefix, bool is_last) {

  // 1. Draw the current node's branch
  std::cout << prefix;
  if (nodes[node].parent_idx != -1) { // If not the root
    std::cout << (is_last ? "└── " : "├── ");
  }

  // 2. Print the node ID (1-indexed for math notation) and its parameters
  std::cout << "x" << (node + 1);
  if (nodes[node].parent_idx == -1) {
    std::cout << " [P(x" << (node + 1) << "=1) = " << nodes[node].p_given_1
              << "]\n";
  } else {
    int parent = nodes[node].parent_idx;
    std::cout << " [P(1|x" << (parent + 1) << "=0)=" << nodes[node].p_given_0
              << ", P(1|x" << (parent + 1) << "=1)=" << nodes[node].p_given_1
              << "]\n";
  }

  // 3. Calculate the prefix for the next depth level
  std::string child_prefix = prefix;
  if (nodes[node].parent_idx != -1) {
    child_prefix += is_last ? "    " : "│   ";
  }

  // 4. Recurse for all children
  for (size_t i = 0; i < children[node].size(); ++i) {
    bool child_is_last = (i == children[node].size() - 1);
    PrintTreeRecursive(children[node][i], children, nodes, child_prefix,
                       child_is_last);
  }
}

// Main visualization function
inline void VisualizeTree(const DirectedTree &tree) {
  std::vector<std::vector<int>> children(WORD_SIZE);
  int root = -1;

  // Build a parent-to-child adjacency list from the flat node array
  for (int i = 0; i < WORD_SIZE; ++i) {
    int parent = tree.nodes[i].parent_idx;
    if (parent == -1) {
      root = i;
    } else {
      children[parent].push_back(i);
    }
  }

  if (root != -1) {
    PrintTreeRecursive(root, children, tree.nodes, "", true);
  }
}

inline void
PrintChainProbabilities(const std::array<size_t, WORD_SIZE> &MIMIC_ORDER,
                        const EntropyMetrics &em) {
  uint8_t root = MIMIC_ORDER[0];
  double p_root = static_cast<double>(em.marginal_counts[root]) / DATABASE_SIZE;

  std::cout << "P(x" << (root + 1) << " = 1) = " << p_root << "\n";

  // Print the Conditional Nodes
  for (int i = 1; i < WORD_SIZE; i++) {
    uint8_t child = MIMIC_ORDER[i];
    uint8_t parent = MIMIC_ORDER[i - 1];

    int offset = ((child * WORD_SIZE) + parent) * 4;

    int count_00 = em.joint_counts[offset + 0];
    int count_01 = em.joint_counts[offset + 1];
    int count_10 = em.joint_counts[offset + 2];
    int count_11 = em.joint_counts[offset + 3];

    int parent_0_total = count_00 + count_10;
    int parent_1_total = count_01 + count_11;

    double p_given_0 = (parent_0_total > 0)
                           ? static_cast<double>(count_10) / parent_0_total
                           : 0.0;

    double p_given_1 = (parent_1_total > 0)
                           ? static_cast<double>(count_11) / parent_1_total
                           : 0.0;

    std::cout << "P(x" << (child + 1) << " = 1 | x" << (parent + 1)
              << " = 0) = " << p_given_0 << "\n";
    std::cout << "P(x" << (child + 1) << " = 1 | x" << (parent + 1)
              << " = 1) = " << p_given_1 << "\n";
  }
}

inline DirectedTree
MimicToTree(const std::array<size_t, WORD_SIZE> &mimic_order,
            const EntropyMetrics &em) {
  DirectedTree result;

  int root = mimic_order[0];
  result.eval_order[0] = root;
  result.nodes[root].parent_idx = -1;
  result.nodes[root].p_given_0 = 0.0;
  result.nodes[root].p_given_1 =
      static_cast<double>(em.marginal_counts[root]) / DATABASE_SIZE;

  for (int i = 1; i < WORD_SIZE; ++i) {
    int child = mimic_order[i];
    int parent = mimic_order[i - 1];

    result.eval_order[i] = child;
    result.nodes[child].parent_idx = parent;

    int offset = ((child * WORD_SIZE) + parent) * 4;
    int count_00 = em.joint_counts[offset + 0];
    int count_01 = em.joint_counts[offset + 1];
    int count_10 = em.joint_counts[offset + 2];
    int count_11 = em.joint_counts[offset + 3];

    int parent_0_total = count_00 + count_10;
    int parent_1_total = count_01 + count_11;

    result.nodes[child].p_given_0 =
        (parent_0_total > 0) ? static_cast<double>(count_10) / parent_0_total
                             : 0.0;
    result.nodes[child].p_given_1 =
        (parent_1_total > 0) ? static_cast<double>(count_11) / parent_1_total
                             : 0.0;
  }

  return result;
}

inline DirectedTree IndependentToTree(const EntropyMetrics& em) {
  DirectedTree result;
  
  for (int i = 0; i < WORD_SIZE; ++i) {
    result.eval_order[i] = i;              // Order doesn't matter for independent variables
    result.nodes[i].parent_idx = -1;       // -1 means no parent (root node)
    result.nodes[i].p_given_0 = 0.0;       // Unused when parent_idx is -1
    result.nodes[i].p_given_1 = static_cast<double>(em.marginal_counts[i]) / DATABASE_SIZE;
  }
  
  return result;
}
