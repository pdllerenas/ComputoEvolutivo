#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../include/evolution/Grid1D.hpp"
#include "../include/evolution/KFLS.hpp"
#include "../include/problem/DIMACInstance.hpp"
#include "../include/problem/FrequencyAssignment.hpp"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string filename = "kfls_results_rank_" + std::to_string(rank) + ".csv";
  std::ofstream outfile(filename);
  outfile << "k,run_id,penalty\n"; // CSV Header
  outfile.close();

  std::string conv_filename =
      "convergence_rank_" + std::to_string(rank) + ".csv";
  std::ofstream conv_file(conv_filename);
  conv_file << "k,run_id,generation,best_penalty\n"; // Header
  conv_file.close();

  if (argc != 4) {
    std::cerr << "Invalid Argument Count.\nUsage: " << argv[0]
              << " <population_size> <max_generations> <dimacs_filepath>\n";
    MPI_Finalize();
    return 1;
  }

  size_t pop_size = std::atoll(argv[1]);
  size_t max_generations = std::atoll(argv[2]);
  std::string instance_file = argv[3];

  const int NUM_RUNS = 10;

  DIMACSInstance graph;

  try {
    graph = ParseDIMACS(instance_file);
  } catch (const std::exception &e) {
    if (rank == 0)
      std::cerr << "Rank " << rank << " Error: " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  size_t d = graph.num_nodes;

  size_t total_k_values = pop_size;
  size_t chunk_size = total_k_values / size;
  size_t remainder = total_k_values % size;
  size_t start_k =
      rank * chunk_size + std::min(static_cast<size_t>(rank), remainder) + 1;
  size_t end_k = start_k + chunk_size + (rank < remainder ? 1 : 0) - 1;

  // Safe upper bound is d (number of nodes), without any pre-computing
  FrequencyAssignmentEvaluator evaluator(0, d, graph.interference_edges);

  std::vector<double> local_k_means(end_k - start_k + 1);
  std::vector<int> local_k_values(end_k - start_k + 1);

  for (size_t k = start_k; k <= end_k; ++k) {
    size_t local_idx = k - start_k;
    local_k_values[local_idx] = k;

    std::vector<double> run_penalties(NUM_RUNS);
    std::vector<std::vector<double>> run_convergences(NUM_RUNS);

#pragma omp parallel for schedule(dynamic)
    for (int run = 0; run < NUM_RUNS; ++run) {
      double run_start_time = MPI_Wtime();

      std::mt19937 rng(42 + run + (k * 1000));
      Grid1D<int> local_grid(pop_size, d);

      local_grid.Initialize(evaluator.GetLowerBound(),
                            evaluator.GetUpperBound(), rng);
      KFLS<int> local_sweeper(local_grid, k, max_generations);

      auto [best_genes, history] =
          local_sweeper.Minimize(evaluator, rng, run_start_time, 300.0);
      run_penalties[run] = evaluator.Evaluate(best_genes);
      run_convergences[run] = std::move(history);
    }

    std::ofstream append_file(filename, std::ios_base::app);
    for (int run = 0; run < NUM_RUNS; ++run) {
      append_file << k << "," << run << "," << run_penalties[run] << "\n";
    }
    append_file.close();

    std::ofstream append_conv(conv_filename, std::ios_base::app);
    for (int run = 0; run < NUM_RUNS; ++run) {
      for (size_t gen = 0; gen < run_convergences[run].size(); ++gen) {
        append_conv << k << "," << run << "," << gen << ","
                    << run_convergences[run][gen] << "\n";
      }
    }
    append_conv.close();

    double sum =
        std::accumulate(run_penalties.begin(), run_penalties.end(), 0.0);
    local_k_means[local_idx] = sum / NUM_RUNS;

    std::cout << "Rank " << rank << " finished k=" << k
              << " | Mean Penalty: " << local_k_means[local_idx] << "\n";
  }

  std::vector<int> recvcounts(size);
  int local_count = local_k_means.size();
  MPI_Gather(&local_count, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  std::vector<int> displs(size, 0);
  std::vector<double> global_k_means;
  std::vector<int> global_k_values;

  if (rank == 0) {
    global_k_means.resize(total_k_values);
    global_k_values.resize(total_k_values);
    for (int i = 1; i < size; ++i) {
      displs[i] = displs[i - 1] + recvcounts[i - 1];
    }
  }

  MPI_Gatherv(local_k_means.data(), local_count, MPI_DOUBLE,
              global_k_means.data(), recvcounts.data(), displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_k_values.data(), local_count, MPI_INT,
              global_k_values.data(), recvcounts.data(), displs.data(), MPI_INT,
              0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "\n--- Final k-Parameter Benchmark ---\n";
    for (size_t i = 0; i < total_k_values; ++i) {
      std::cout << "k = " << std::setw(3) << global_k_values[i]
                << " | Mean Penalty: " << std::fixed << std::setprecision(4)
                << global_k_means[i] << "\n";
    }
  }

  MPI_Finalize();

  return 0;
}
