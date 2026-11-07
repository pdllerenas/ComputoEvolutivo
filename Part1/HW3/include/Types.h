#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

/*
 * @brief Node dependency and probability values
 *
 */
struct DependencyNode {
  int16_t parent_idx; /**< Index value of parent. -1 for root */
  double p_given_0;   /**< Probability of P(X|X_parent = 0). 0 for root */
  double p_given_1; /**< Probability of P(X|X_parent = 1). Marginal for root */
};

/*
 * Dependency nodes for DB1:
 * x1 -> x2 -> x3 -> x4 -> x5 -> x6
 */
constexpr std::array<DependencyNode, 6> CHAIN = {{
    {-1, 0.20, 0.00}, // x1 Root
    {0, 0.92, 0.31},  // x2 Parent is x1
    {1, 0.42, 0.15},  // x3 Parent is x2
    {2, 0.28, 0.75},  // x4 Parent is x3
    {3, 0.67, 0.43},  // x5 Parent is x4
    {4, 0.35, 0.82}   // x6 Parent is x5
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
    {5, 0.42, 0.15},  // x3 Parent is x6
    {0, 0.28, 0.75},  // x4 Parent is x1
    {1, 0.67, 0.43},  // x5 Parent is x2
    {1, 0.35, 0.82}   // x6 Parent is x2
}};

// Topological order to evaluate TREE
constexpr std::array<uint16_t, 6> TREE_EVAL_ORDER = {0, 1, 3, 4, 5, 2};

/*
 * Database generation methods
 */
constexpr size_t DATABASE_SIZE = 100'000;
constexpr size_t WORD_SIZE = 6;

struct EntropyMetrics {
  std::vector<int> marginal_counts;
  std::vector<int> joint_counts;
};

struct DisjointSets {
  std::vector<int> parent, rnk;
  int n;

  DisjointSets(int _n) : n(_n), parent(_n + 1), rnk(_n + 1) {
    for (int i = 0; i <= n; ++i) {
      rnk[i] = 0;
      parent[i] = i;
    }
  }

  int find(int u) {
    if (u != parent[u]) {
      parent[u] = find(parent[u]);
    }
    return parent[u];
  }

  void merge(int x, int y) {
    x = find(x), y = find(y);
    if (rnk[x] > rnk[y])
      parent[y] = x;
    else
      parent[x] = y;

    if (rnk[x] == rnk[y])
      rnk[y]++;
  }
};

struct Graph {
  int V, E;
  std::vector<std::pair<double, std::pair<int, int>>> edges;
  Graph(int _V, int _E) : V(_V), E(_E) {}

  void addEdge(int u, int v, double w) { edges.push_back({w, {u, v}}); }

  /*
   * Finds the Maximum Spanning Tree
   */
  int kruskalMST() {
    int mst_wt = 0;
    // sort in descending order
    std::sort(edges.rbegin(), edges.rend());
    DisjointSets ds(V);

    for (auto it = edges.begin(); it != edges.end(); ++it) {
      int u = it->second.first;
      int v = it->second.second;

      int set_u = ds.find(u);
      int set_v = ds.find(v);

      if (set_u != set_v) {
        mst_wt += it->first;
        ds.merge(set_u, set_v);
      }
    }
    return mst_wt;
  }

  std::vector<std::pair<int, int>> extractMST() {
    std::vector<std::pair<int, int>> mst_edges;
    std::sort(edges.rbegin(), edges.rend());
    DisjointSets ds(V);

    for (const auto &edge : edges) {
      int u = edge.second.first;
      int v = edge.second.second;

      int set_u = ds.find(u);
      int set_v = ds.find(v);

      if (set_u != set_v) {
        mst_edges.push_back({u, v});
        ds.merge(set_u, set_v);
      }
    }
    return mst_edges;
  }

  void print() const noexcept {
    std::cout << "Graph (V=" << V << ", E=" << E << ")\n";
    for (const auto &[weight, nodes] : edges) {
      const auto &[u, v] = nodes;
      std::cout << (u + 1) << " -- " << (v + 1) << " [weight: " << weight
                << "]\n";
    }
  }
};

struct DirectedTree {
  std::array<DependencyNode, WORD_SIZE> nodes;
  std::array<uint16_t, WORD_SIZE> eval_order;
};
