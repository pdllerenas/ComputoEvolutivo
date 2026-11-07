#pragma once

#include <cmath>
#include <vector>

#include "Types.h"

struct EpsilonCompare {
  const std::vector<double> &fitness;
  const std::vector<double> &violation;
  double eps = 1e-6;

  EpsilonCompare(const std::vector<double> &fit,
                 const std::vector<double> &viol, double e)
      : fitness(fit), violation(viol), eps(e) {}

  bool operator()(size_t idx1, size_t idx2) const {
    double f1 = fitness[idx1];
    double f2 = fitness[idx2];
    double v1 = violation[idx1];
    double v2 = violation[idx2];

    if (v1 <= eps && v2 <= eps) {
      return f1 < f2;
    } else if (v1 <= eps) {
      return true;
    } else if (v2 <= eps) {
      return false;
    } else {
      return v1 < v2;
    }
  }
};

class Takahama {
private:
  double eps0;
  double cp;
  size_t T_control;

public:
  Takahama(double initial_eps, size_t control_gens, double decay_cp = 2.0)
      : eps0(initial_eps), T_control(control_gens), cp(decay_cp) {}
  double operator()(size_t current_gen) const {
    if (current_gen >= T_control) {
      return 0.0;
    }
    double fraction = 1.0 - (static_cast<double>(current_gen) /
                             static_cast<double>(T_control));

    return eps0 * std::pow(fraction, cp);
  }
};

static void batch_evaluate(Population &pop, const Objective &obj,
                    const HimmelblauConstraints &constraints) {
  pop.all_fitnesses.resize(pop.num_individuals);
  pop.all_violations.resize(pop.num_individuals);

// #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < pop.num_individuals; ++i) {
    const double *individual_genes = &pop.all_genes[i * pop.gene_size];

    pop.all_fitnesses[i] = obj(individual_genes);
    pop.all_violations[i] = constraints(individual_genes);
  }
}

static double get_initial_eps0(const std::vector<double> &all_violations) {
  double max_v = 0.0;

  // #pragma omp parallel for reduction(max:max_v)
  // for (size_t i = 0; i < all_violations.size(); ++i) {
  //     if (all_violations[i] > max_v) {
  //         max_v = all_violations[i];
  //     }
  // }

  max_v = *std::max_element(all_violations.begin(), all_violations.end());

  return max_v;
}
