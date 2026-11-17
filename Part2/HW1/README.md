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
mpirun -np 4 ./src/main <population_size> <max_generations> <data.csv>
```

### Example Run

To solve the Frequency Assignment problem with a population of 50 for 1,000 generations, for some Graph Coloring Problem:

```bash
./src/main 50 1000 file.csv
```

### Output Interpretation

The program will execute 10 independent runs (with different seeds). It will output the best penalty
    found in each run, followed by the global aggregated statistics.
