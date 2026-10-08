#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

struct Population {
  std::vector<double> all_genes;
  std::vector<double> all_fitnesses;
  size_t num_individuals = 0;
  size_t gene_size = 0;
};

/*
 * @brief Node dependency and probability values
 *
 */
struct DependencyNode {
  int16_t parent_idx; /**< Index value of parent. -1 for root */
  double p_given_0;   /**< Probability of P(X|X_parent = 0). 0 for root */
  double p_given_1; /**< Probability of P(X|X_parent = 1). Marginal for root */
};

constexpr size_t DATABASE_SIZE = 1000;

template <typename T> struct DisjointSets {
  std::vector<T> parent, rnk;
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

  void merge(T x, T y) {
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
    DisjointSets<int> ds(V);

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
    DisjointSets<int> ds(V);

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
  std::vector<DependencyNode> nodes;
  std::vector<uint16_t> eval_order;
};

class Objective {
public:
  virtual double operator()(const std::vector<double> &x) const = 0;
  virtual std::pair<double, double> bounds() const = 0;
  virtual size_t dimension() const = 0;
  virtual double optimal() const = 0;
  virtual ~Objective() = default;
};

// Global minimum for Ackley: f(x) = 0 at x = (0,...,0)
class Ackley : public Objective {
private:
  double a, b, c;
  size_t dim;

public:
  Ackley(size_t dim_ = 2, double a_ = 20.0, double b_ = 0.2,
         double c_ = 2 * std::numbers::pi)
      : dim(dim_), a(a_), b(b_), c(c_) {}

  double operator()(const std::vector<double> &x) const override {
    int d = x.size();
    double sum_sq = 0.0, sum_cos = 0.0;
    for (double y : x) {
      sum_sq += y * y;
      sum_cos += std::cos(c * y);
    }
    return -a * std::exp(-b * std::sqrt(sum_sq / d)) - std::exp(sum_cos / d) +
           a + std::numbers::e;
  }

  size_t dimension() const { return dim; }

  std::pair<double, double> bounds() const {
    return std::make_pair<double, double>(-32.768, 32.768);
  }

  double optimal() const { return 0.0; }

  ~Ackley() = default;
};

// Global minimum for Ackley: f(x) = 0 at x = (0,...,0)
class Griewank : public Objective {
private:
  size_t dim;

public:
  Griewank(size_t dim_ = 2) : dim(dim_) {}
  double operator()(const std::vector<double> &x) const override {
    int d = x.size();
    double sum_sq = 0.0, prod_cos = 1.0;
    for (size_t i = 0; i < d; ++i) {
      double y = x[i];
      sum_sq += y * y;
      prod_cos *= std::cos(y / std::sqrt(i + 1));
    }
    return sum_sq / 4000.0 - prod_cos + 1.0;
  }

  size_t dimension() const { return dim; }

  std::pair<double, double> bounds() const {
    return std::make_pair<double, double>(-600.0, 600.0);
  }

  double optimal() const { return 0.0; }

  ~Griewank() = default;
};

// Global minimum for Sphere: f(x) = 0 at x = (0,...,0)
class Sphere : public Objective {
private:
  size_t dim;

public:
  Sphere(size_t dim_ = 2) : dim(dim_) {}
  double operator()(const std::vector<double> &x) const override {
    int d = x.size();
    double sum_sq = 0.0;
    for (double y : x) {
      sum_sq += y * y;
    }
    return sum_sq;
  }

  size_t dimension() const { return dim; }

  std::pair<double, double> bounds() const {
    return std::make_pair<double, double>(-5.12, 5.12);
  }

  double optimal() const { return 0.0; }

  ~Sphere() = default;
};

struct EDAResult {
	std::vector<double> best_gene;
	double best_fitness;
	size_t generations_used;
	bool converged;
};
