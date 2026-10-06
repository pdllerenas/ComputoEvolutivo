#include <iostream>
#include <random>

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
    unsigned semilla = rd();
    std::mt19937 gen(semilla);

    Population pop = RandomPopulation(obj, n, gen);
    std::vector<double> resultado = EMNA(pop, obj, gen, MAX_EVALS);
    double mejor_valor = obj(resultado);

    std::cerr << "  run " << (i + 1) << "/" << num_runs
              << " (seed=" << semilla << ") -> " << mejor_valor << "\n";

    best_vals.push_back(mejor_valor);
  }
  return best_vals;
}

int main(int argc, char **argv) {
  std::random_device rd;
  std::mt19937 rng(rd());

  constexpr size_t DIM = 10;

  std::vector<double> means(10, 0.0);
  std::vector<double> devs(10, 1.0);

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

  }

}
