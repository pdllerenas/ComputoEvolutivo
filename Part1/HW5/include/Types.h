#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

struct Population {
  std::vector<double> all_genes;
  std::vector<double> all_fitnesses;
  std::vector<double> all_violations;
  size_t num_individuals = 0;
  size_t gene_size = 0;
};

class Objective {
public:
  virtual double operator()(const double *x) const = 0;
  virtual std::pair<double, double> range() const = 0;
  virtual size_t dimension() const = 0;
  virtual double optimal() const = 0;
  virtual ~Objective() = default;
};

class Himmelblau : public Objective {
private:
  size_t dim;

public:
  Himmelblau(size_t dim_) : dim(dim_) {}

  double operator()(const double *x) const override {
    return (x[0] * x[0] + x[1] - 11) * (x[0] * x[0] + x[1] - 11) +
           (x[0] + x[1] * x[1] - 7) * (x[0] + x[1] * x[1] - 7);
  }
  size_t dimension() const override { return dim; }
  std::pair<double, double> range() const override {
    return std::make_pair(0.0, HUGE_VAL);
  }
  double optimal() const override { return 0.0; }
  ~Himmelblau() = default;
};

struct HimmelblauConstraints {
  double operator()(const double *x) const {
    double g1 = (x[0] - 5.0) * (x[0] - 5.0) + x[1] * x[1] - 26.0;
    double g2 = 4.0 * x[0] + x[1] - 20.0;
    return std::max(0.0, g1) + std::max(0.0, g2);
  }

  std::pair<double, double> inspect(const std::vector<double> &x) const {
    double g1 = (x[0] - 5.0) * (x[0] - 5.0) - x[1] * x[1] - 26.0;
    double g2 = 4.0 * x[0] + x[1] - 20.0;
    return {g1, g2};
  }
};
