#pragma once

#include <cmath>
#include <utility>
#include <vector>

class FrequencyAssignmentEvaluator {
private:
  int lower_bound;
  int upper_bound;
  std::vector<std::pair<int, int>> interference_edges;

public:
  FrequencyAssignmentEvaluator(int min_freq, int max_freq,
                               const std::vector<std::pair<int, int>> &edges)
      : lower_bound(min_freq), upper_bound(max_freq),
        interference_edges(edges) {}

  [[nodiscard]] int GetLowerBound() const { return lower_bound; }
  [[nodiscard]] int GetUpperBound() const { return upper_bound; }

  [[nodiscard]] double Evaluate(const std::vector<int> &frequencies) const {
    double total_penalty = 0.0;

    for (const auto &edge : interference_edges) {
      int node_u = edge.first;
      int node_v = edge.second;

			int freq_diff = std::abs(frequencies[node_u] - frequencies[node_v]);

			if (freq_diff == 0) {
				total_penalty += 100.0;
			} else if (freq_diff == 1) {
				total_penalty += 10.0;
			}
    }
		return total_penalty;
  }
};
