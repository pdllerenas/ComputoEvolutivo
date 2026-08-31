# Computo Evolutivo

## EMNA
To run EMNA, you must have `cmake` and `make`. Run
```sh
mkdir build && cd build && cmake .. && make
```
This will create a binary file which can be executed with
```sh
./src/emna <population size> <dimension>
```

For example,
```sh
./src/emna 1000 2
```
runs on 1000 population size, on three benchmark functions at dimension 2.

To prevent spam on the terminal, re-direct the (stderr) output to a file:
```sh
./src/emna 1000 2 2> experiment.log
```
