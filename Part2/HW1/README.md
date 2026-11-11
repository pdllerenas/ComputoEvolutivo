# High-Performance Cellular Genetic Algorithm (kFLS)

A high-performance, cache-friendly implementation of a Cellular Genetic
Algorithm utilizing the **k-Fixed Line Sweep (kFLS)** update policy. Written in
standard **C++20**, this project is designed to solve discrete, multi-modal
optimization problems—including the NP-hard Radio Frequency Assignment Problem
(FAP).

## Features

* **k-Fixed Line Sweep (kFLS) Policy**: Bridges the gap between synchronous and
  asynchronous cellular updates. Control the propagation of genetic material by
  tuning the $k$-block size.
* **Cache-Optimal Memory Architecture**: Uses a flat, contiguous 1D
  `std::vector` layout with double-buffering (read/write buffers) managed by a
  central `Grid1D` topology.
* **Zero Inner-Loop Allocations**: Operators are decoupled and iterator-based.
  Crossover writes directly into pre-allocated memory buffers, completely
  eliminating heap allocations (`new`/`malloc`) during the evolutionary sweep.
* **Plug-and-Play Evaluators**: Strictly separates the evolutionary engine from
  the problem domain. Swap between continuous benchmark functions (like
  Discrete Ackley) and NP-hard graph problems (like Frequency Assignment)
  without touching the algorithm's core.
* **Statistical Benchmarking**: Automatically wraps execution in a 30-run
  statistical aggregator using isolated, reproducible random seeds to evaluate
  algorithmic robustness.

## Project Structure

```text
.
├── CMakeLists.txt                # C++20, -O3, -march=native build configuration
├── include/
│   ├── evolution/
│   │   ├── Grid1D.hpp            # Spatial topology and double-buffer management
│   │   ├── KFLS.hpp              # The overarching k-Fixed Line Sweep orchestrator
│   │   └── Population.hpp        # Contiguous 1D memory store
│   ├── operators/
│   │   ├── Crossover.hpp         # Iterator-based Bounded Discrete Crossover
│   │   ├── Mutation.hpp          # In-place Uniform Mutation
│   │   └── Selection.hpp         # L-Best East-West Tournament Selection with tie-breaking
│   └── problem/
│       └── FrequencyAssignment.hpp # NP-hard benchmark evaluator
└── src/
    └── main.cpp                  # CLI parser, dependency injector, and statistical runner
```

## Requirements

* **C++20 Compatible Compiler**: GCC 10+, Clang 11+, or MSVC 19.28+
* **CMake**: Version 3.14 or higher

## Build Instructions

This project is configured for aggressive compiler optimizations (`-O3`, `-march=native`, `-flto`). To build the executable:

```bash
# 1. Create a build directory
mkdir build && cd build

# 2. Generate build files
cmake ..

# 3. Compile the project
make
```

## Usage

The executable takes four positional arguments:
```bash
./src/main <population_size> <dimension> <k_block_size> <max_generations>
```

* **`population_size`**: Total number of individuals in the 1D grid.
* **`dimension`**: Number of genes per individual (e.g., number of nodes in the Frequency Assignment graph).
* **`k`**: The block size for the kFLS sweep. 
  * Set `k = 1` for purely asynchronous updates.
  * Set `k = population_size` for purely synchronous updates.
  * Set `1 < k < population_size` to balance exploitation and exploration.
* **`max_generations`**: The stopping criterion.

### Example Run

To solve a 50-node Frequency Assignment problem with a population of 100,
updating in k-blocks of 10 for 1,000 generations:

```bash
./src/main 100 50 10 1000
```

### Output Interpretation

The program will execute 30 independent runs (seeds 42 through 71) to account
for the stochastic nature of metaheuristics. It will output the best penalty
    found in each run, followed by the global aggregated statistics:

```text
Starting k-Fixed Line Sweep (kFLS) Cellular GA...
Population: 100 | k-Block: 10 | Runs: 30

Run  1 | Best Penalty: 0
Run  2 | Best Penalty: 10
...
Run 30 | Best Penalty: 0

--- Final Statistics (30 Runs) ---
Global Best : 0
Mean Penalty: 2.3333
Std Dev     : 4.3018
```

