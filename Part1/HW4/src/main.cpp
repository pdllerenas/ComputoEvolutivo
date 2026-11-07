#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

#include "../include/ChowLiu.h"
#include "../include/MIMIC.h"
#include "../include/Types.h"

double ComputeStdDev(const std::vector<double> &data, double mean) {
  if (data.size() < 2)
    return 0.0;

  double sum_of_squares = std::accumulate(data.begin(), data.end(), 0.0,
                                          [mean](double acc, double val) {
                                            double diff = val - mean;
                                            return acc + (diff * diff);
                                          });

  return std::sqrt(sum_of_squares / (data.size() - 1));
}

struct Entry {
  std::string name;
  const Objective *obj;
  size_t max_evals;
};

template <typename func>
std::tuple<std::vector<double>, size_t>
Experiment(const Objective &obj, func EDA, size_t n, size_t MAX_EVALS,
           size_t num_runs = 10) {
  std::vector<double> best_vals;
  best_vals.reserve(num_runs);
  size_t mean_iterations = 0;

  std::random_device rd;

  for (int i = 0; i < num_runs; ++i) {
    unsigned seed = rd();
    std::mt19937 rng(seed);

    EDAResult result = EDA(obj, n, MAX_EVALS, rng);
    best_vals.push_back(result.best_fitness);
    mean_iterations += result.generations_used;
  }
  return std::make_tuple(best_vals,
                         static_cast<double>(mean_iterations) / num_runs);
}

template <typename func>
void RunAndReport(const std::string &algo_name, const Objective &obj, func EDA,
                  size_t n, size_t max_evals, int num_runs = 10) {

  auto [best_vals, mean_iters] = Experiment(obj, EDA, n, max_evals, num_runs);

  std::sort(best_vals.begin(), best_vals.end());

  size_t count = best_vals.size();
  double best = best_vals.front();
  double worst = best_vals.back();

  double mean =
      std::accumulate(best_vals.begin(), best_vals.end(), 0.0) / count;

  double median = (count % 2 == 0)
                      ? (best_vals[count / 2 - 1] + best_vals[count / 2]) / 2.0
                      : best_vals[count / 2];

  double stddev = ComputeStdDev(best_vals, mean);

  std::cout << algo_name << "," << best << "," << worst << "," << mean << ","
            << median << "," << stddev << "," << mean_iters << "\n";
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "Invalid arguments. Usage: " << argv[0]
              << " <population> <dimension>\n";
    return 1;
  }

  size_t n = static_cast<size_t>(std::atoi(argv[1]));
  size_t d = static_cast<size_t>(std::atoi(argv[2]));

  Ackley ackley(d);
  Griewank griewank(d);
  Sphere sphere(d);

  std::vector<Entry> problems = {
      {"Ackley", &ackley, 10000 * ackley.dimension()},
      {"Griewank", &griewank, 10000 * griewank.dimension()},
      {"Sphere", &sphere, 10000 * sphere.dimension()},
  };

  for (const auto &problem : problems) {
    std::cout << "\n=== " << problem.name << " ===\n";
    std::cout << "algorithm,best,worst,mean,median,stddev,mean_generations\n";

    RunAndReport("Chow Liu", *problem.obj, ChowLiu::EDA, n, problem.max_evals,
                 10);
    RunAndReport("MIMIC", *problem.obj, MIMIC::EDA, n, problem.max_evals, 10);
  }

  return 0;
}
