#include <iostream>
#include <random>

#include "../include/ChowLiu.h"
#include "../include/Generator.h"
#include "../include/Selection.h"
#include "../include/Stats.h"
#include "../include/Types.h"

struct Entry {
  std::string name;
  const Objective *obj;
  size_t max_evals;
};

std::vector<double> Experiment(const Objective &obj, size_t n, size_t MAX_EVALS,
                               int num_runs = 10) {
  std::vector<double> best_vals;
  best_vals.reserve(num_runs);

  std::random_device rd;

  for (int i = 0; i < num_runs; ++i) {
    unsigned seed = rd();
    std::mt19937 rng(seed);

    EDAResult result = ChowLiu::EDA(obj, n, MAX_EVALS, rng);
    best_vals.push_back(result.best_fitness);
  }
  return best_vals;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "Invalid arguments. Usage: " << argv[0]
              << " <population> <dimension>\n";
  }

  size_t n = static_cast<size_t>(std::atoi(argv[1]));
  size_t d = static_cast<size_t>(std::atoi(argv[2]));

  Ackley ackley(d);
  Griewank griewank(d);
  Sphere sphere(d);

  std::cout << "run,problem,final_error" << std::endl;

  std::vector<Entry> problems = {
      {"Ackley", &ackley, 100000 * ackley.dimension()},
      {"Griewank", &griewank, 100000 * griewank.dimension()},
      {"Sphere", &sphere, 100000 * sphere.dimension()},
  };

  for (const auto &problem : problems) {
    std::vector<double> best_vals =
        Experiment(*problem.obj, n, problem.max_evals, 10);
    for (size_t i = 0; i < best_vals.size(); ++i) {
      std::cout << i << ", " << problem.name << ", "
                << best_vals[i] - std::abs((*problem.obj).optimal()) << '\n';
    }
  }
}
