#pragma once

#include <algorithm>
#include <utility>
#include <limits>
#include <mpi.h>
#include <random>
#include <vector>

#include "../operators/Crossover.hpp"
#include "../operators/Mutation.hpp"
#include "../operators/Selection.hpp"
#include "Grid1D.hpp"

template <typename T> class KFLS {
private:
  Grid1D<T> grid;
  size_t k;
  size_t max_generations;

public:
  KFLS(Grid1D<T> &_grid, size_t block_size, size_t _max_generations)
      : grid(_grid), k(block_size), max_generations(_max_generations) {}

  // Evaluator must provide the following:
  // double Evaluate(const std::vector<T> &genes)
  // GetLowerBound()
  // GetUpperBound()
  template <typename Evaluator>
  std::pair<std::vector<T>, std::vector<double>> Minimize(const Evaluator &evaluator, std::mt19937 &rng,
                          double global_start_time, double timeout_seconds) {
    size_t pop_size = grid.GetPopulationSize();
    size_t d = grid.GetGeneSize();

    T lower_bound = evaluator.GetLowerBound();
    T upper_bound = evaluator.GetUpperBound();

    std::vector<T> best_individual(d);
    double best_fitness = std::numeric_limits<double>::infinity();

    std::vector<double> convergence_history;
    convergence_history.reserve(max_generations);

    for (size_t i = 0; i < pop_size; ++i) {
      std::vector<T> current_genes(grid.GetReadGenesBegin(i),
                                   grid.GetReadGenesBegin(i) + d);
      double fit = evaluator.Evaluate(current_genes);
      grid.SetReadFitness(i, fit);

      bool is_better = (fit < best_fitness);
      if (is_better) {
        best_fitness = fit;
        best_individual = current_genes;
      }
    }

    for (size_t gen = 0; gen < max_generations; ++gen) {
      for (size_t start_idx = 0; start_idx < pop_size; start_idx += k) {
        size_t block_end = std::min(start_idx + k, pop_size);
        for (size_t i = start_idx; i < block_end; ++i) {
          size_t west_idx = grid.GetWestNeighbor(i);
          size_t east_idx = grid.GetEastNeighbor(i);

          size_t partner =
              EastWestSelection<T>(i, west_idx, east_idx, grid, rng);

          auto mother_iterator = grid.GetReadGenesBegin(i);
          auto father_iterator = grid.GetReadGenesBegin(partner);

          std::vector<T> child = BoundedDiscreteCrossover(
              mother_iterator, father_iterator, d, rng);
					double mutation_rate = 1.0 / d;
          UniformMutation(child, mutation_rate, lower_bound, upper_bound, rng);

          double child_fitness = evaluator.Evaluate(child);
          double mother_fitness = grid.GetReadFitness(i);

          bool replace = (child_fitness < mother_fitness);

          if (replace) {
            grid.SetWriteIndividual(i, child, child_fitness);
            bool is_global_best = (child_fitness < best_fitness);

            if (is_global_best) {
              best_fitness = child_fitness;
              best_individual = child;
            }
          } else {
            grid.SetWriteIndividual(i, mother_iterator, mother_fitness);
          }
        }
        grid.FlushKBlocks(start_idx, k);
      }
      convergence_history.push_back(best_fitness);
      if (gen % 20 == 0) {
        if (MPI_Wtime() - global_start_time > timeout_seconds) {
          break;
        }
      }
    }
    return {best_individual, convergence_history};
  }
};
