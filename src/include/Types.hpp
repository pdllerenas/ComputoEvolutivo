#pragma once

#include <cmath>
#include <compare>
#include <numbers>
#include <vector>

struct Individual {
  std::vector<double> genes;
  double fitness = 0.0;
  size_t dim;

  Individual(const std::vector<double> &_genes, double _fitness)
      : genes(_genes), fitness(_fitness) {
    dim = _genes.size();
  }

  std::partial_ordering operator<=>(const Individual &other) const {
    return fitness <=> other.fitness;
  }

  std::vector<double> GetGenes() const { return genes; }
  size_t GetDimension() const { return dim; }

  friend void swap(Individual &lhs, Individual &rhs) noexcept {
    std::swap(lhs.genes, rhs.genes);
    std::swap(lhs.fitness, rhs.fitness);
  }
};

// This arrangement allows for cache-friendly reads, as genes are contiguous
struct Population {
  std::vector<double> all_genes;
  std::vector<double> all_fitnesses;
  size_t num_individuals = 0;
  size_t gene_size = 0;
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
