#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include "../include/ChowLiu.h"
#include "../include/Stats.h"
#include "../include/Types.h"


Graph BuildChowLiuGraph(const std::array<uint8_t, DATABASE_SIZE> &db,
                 const EntropyMetrics &em) {

  Graph g(WORD_SIZE, WORD_SIZE * (WORD_SIZE - 1) / 2);
  for (uint8_t i = 0; i < WORD_SIZE; ++i) {
    for (uint8_t j = i + 1; j < WORD_SIZE; ++j) {
      int offset = ((i * WORD_SIZE) + j) * 4;
      // I(X,Y) = H(X) - H(X|Y)
      double mutual_info = MarginalEntropy(em.marginal_counts[i]) -
                           ConditionalEntropy(&em.joint_counts[offset]);
      g.addEdge(i, j, mutual_info);
    }
  }
  return g;
}

DirectedTree BuildDirectedTree(Graph &g, const EntropyMetrics &em,
                               std::mt19937 &rng) {
  DirectedTree result;
  auto mst_edges = g.extractMST();

  std::vector<std::vector<int>> adj(WORD_SIZE);
  for (const auto &[u, v] : mst_edges) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  std::uniform_int_distribution<int> dist(0, WORD_SIZE - 1);
  int root = dist(rng);

  bool visited[WORD_SIZE] = {false};
  uint16_t queue[WORD_SIZE];
  int head = 0, tail = 0;

  queue[tail++] = root;
  visited[root] = true;

  result.nodes[root].parent_idx = -1;
  result.nodes[root].p_given_0 = 0.0;
  result.nodes[root].p_given_1 =
      static_cast<double>(em.marginal_counts[root]) / DATABASE_SIZE;

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

        int child = neighbor;
        int parent = current;
        int offset = ((child * WORD_SIZE) + parent) * 4;

        int count_00 = em.joint_counts[offset + 0];
        int count_01 = em.joint_counts[offset + 1];
        int count_10 = em.joint_counts[offset + 2];
        int count_11 = em.joint_counts[offset + 3];

        int parent_0_total = count_00 + count_10;
        int parent_1_total = count_01 + count_11;

        result.nodes[child].p_given_0 =
            (parent_0_total > 0)
                ? static_cast<double>(count_10) / parent_0_total
                : 0.0;
        result.nodes[child].p_given_1 =
            (parent_1_total > 0)
                ? static_cast<double>(count_11) / parent_1_total
                : 0.0;
      }
    }
  }
	return result;
}
