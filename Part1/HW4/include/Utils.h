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

inline void PrintChain(const std::vector<int> &MIMIC_ORDER) {
  for (int i = 0; i < MIMIC_ORDER.size(); ++i) {
    std::cout << (MIMIC_ORDER[i] + 1);
    if (i < MIMIC_ORDER.size() - 1) {
      std::cout << " -> ";
    }
  }
  std::cout << '\n';
}
// Forward declaration for the recursive helper
static void
PrintTreeRecursive(int node, const std::vector<std::vector<int>> &children,
                   const std::vector<DependencyNode> &nodes,
                   const std::vector<double> &means,
                   const std::vector<double> &variances,
                   const std::vector<std::vector<double>> &correlations,
                   std::string prefix, bool is_last) {

  // 1. Draw the current node's branch
  std::cout << prefix;
  if (nodes[node].parent_idx != -1) { // If not the root
    std::cout << (is_last ? "└── " : "├── ");
  }

  // 2. Print the node ID (1-indexed for math notation) and its parameters
  double mu = means[node];
  double sigma = std::sqrt(std::max(0.0, variances[node]));

  std::cout << "x" << (node + 1);
  if (nodes[node].parent_idx == -1) {
    std::cout << " [μ=" << mu << ", σ=" << sigma << "]\n";
  } else {
    int parent = nodes[node].parent_idx;
    double rho = correlations[node][parent];
    std::cout << " [μ=" << mu << ", σ=" << sigma << ", ρ(x" << (parent + 1)
              << ")=" << rho << "]\n";
  }

  // 3. Calculate the prefix for the next depth level
  std::string child_prefix = prefix;
  if (nodes[node].parent_idx != -1) {
    child_prefix += is_last ? "    " : "│   ";
  }

  // 4. Recurse for all children
  for (size_t i = 0; i < children[node].size(); ++i) {
    bool child_is_last = (i == children[node].size() - 1);
    PrintTreeRecursive(children[node][i], children, nodes, means, variances,
                       correlations, child_prefix, child_is_last);
  }
}

// Main visualization function
inline void
VisualizeTree(const DirectedTree &tree, const std::vector<double> &means,
              const std::vector<double> &variances,
              const std::vector<std::vector<double>> &correlations) {
  size_t d = means.size();
  std::vector<std::vector<int>> children(d);
  int root = -1;

  // Build a parent-to-child adjacency list from the flat node array
  for (int i = 0; i < d; ++i) {
    int parent = tree.nodes[i].parent_idx;
    if (parent == -1) {
      root = i;
    } else {
      children[parent].push_back(i);
    }
  }

  if (root != -1) {
    PrintTreeRecursive(root, children, tree.nodes, means, variances,
                       correlations, "", true);
  }
}
