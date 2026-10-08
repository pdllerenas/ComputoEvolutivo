#include <cassert>
#include <random>
#include <vector>

#include "../include/ChowLiu.h"
#include "../include/Generator.h"
#include "../include/Selection.h"
#include "../include/Stats.h"
#include "../include/Types.h"

namespace ChowLiu {
Graph BuildChowLiuGraph(const std::vector<std::vector<double>> &correlations) {
  size_t d = correlations[0].size();
  Graph g(d, d * (d - 1) / 2);
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = i + 1; j < d; ++j) {
      int offset = ((i * d) + j) * 4;
      // I(X,Y) = -1/2 log (1-rho^2)
      double mutual_info = std::abs(correlations[i][j]);
      g.addEdge(i, j, mutual_info);
    }
  }
  return g;
}

DirectedTree BuildDirectedTree(Graph &g, std::mt19937 &rng) {
  DirectedTree result;
  auto mst_edges = g.extractMST();
  size_t d = g.V;

  result.nodes.resize(d);
  result.eval_order.resize(d);

  std::vector<std::vector<int>> adj(d);
  for (const auto &[u, v] : mst_edges) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  std::uniform_int_distribution<int> dist(0, d - 1);
  int root = dist(rng);

  std::vector<bool> visited(d, false);
  std::vector<int> queue(d, false);
  int head = 0, tail = 0;

  queue[tail++] = root;
  visited[root] = true;

  result.nodes[root].parent_idx = -1;

  // bfs
  while (head < tail) {
    int current = queue[head];
    result.eval_order[head] = current;
    head++;

    for (int neighbor : adj[current]) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        queue[tail++] = neighbor;
        result.nodes[neighbor].parent_idx = current;
      }
    }
  }
  return result;
}

EDAResult EDA(const Objective &obj, size_t population_size,
              size_t max_generations, std::mt19937 &rng) {

  Population population = RandomPopulation(obj, population_size, rng);
  size_t n = population.num_individuals;
  size_t d = population.gene_size;
  assert(obj.dimension() == d);
  size_t gen = 0;
  bool converged = false;

  std::vector<double> means;
  std::vector<double> variances;
  std::vector<std::vector<double>> correlations;

  Population elite;

  for (; gen < max_generations; ++gen) {
    elite = SelectElite(population, population_size / 2);

    ComputeMeanVarianceCorrelations(elite, means, variances, correlations);
    // check if the maximum variance across all dimensions is close to 0
    double max_var = *std::max_element(variances.begin(), variances.end());
    if (max_var < 1e-10) {
      converged = true;
      break;
    }
    Graph g = BuildChowLiuGraph(correlations);
    DirectedTree dt = BuildDirectedTree(g, rng);

    // this function replaces population genes and fitnesses, keeping its other
    // attributes intact. Thus, the population variable comes back as
    // the new generation, sampled from the chow-liu tree
    GenerateChowLiuPopulation(population, elite, obj, dt, means, variances,
                              correlations, population_size, rng);
  }
  return {
      std::vector<double>(elite.all_genes.begin(), elite.all_genes.begin() + d),
      elite.all_fitnesses[0], gen, converged};
}

} // namespace ChowLiu
