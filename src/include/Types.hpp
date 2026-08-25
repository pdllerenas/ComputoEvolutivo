#include <vector>

struct Individual {
  std::vector<double> genes;
  double fitness = 0.0;
  std::partial_ordering operator<=>(const Individual &other) const {
    return fitness <=> other.fitness;
  }
};
