#pragma once

#include "Population.hpp"

template <typename T> class Grid1D {
private:
  Population<T> read_buffer;
  Population<T> write_buffer;
  size_t pop_size;
  size_t gene_size;

public:
  Grid1D(size_t size, size_t _gene_size)
      : read_buffer(size, _gene_size), write_buffer(size, _gene_size),
        pop_size(size), gene_size(_gene_size) {}

  void Initialize(T lower, T upper, std::mt19937 &rng) {
    read_buffer.IntializeRandom(lower, upper, rng);
  }

  [[nodiscard]]
  size_t GetWestNeighbor(size_t center_idx) const {
    return (center_idx + pop_size - 1) % pop_size;
  }

  [[nodiscard]]
  size_t GetEastNeighbor(size_t center_idx) const {
    return (center_idx + 1) % pop_size;
  }

  [[nodiscard]]
  auto GetReadGenesBegin(size_t center_idx) const {
    return read_buffer.GetGenesBegin(center_idx);
  }

  void SetReadFitness(size_t target, double fitness) {
    read_buffer.SetFitness(target, fitness);
  }

  [[nodiscard]]
  double GetReadFitness(size_t target) const {
    return read_buffer.GetFitness(target);
  }

  void SetWriteIndividual(size_t idx, const std::vector<T> &genes,
                          double fitness) {
    auto dest_begin = write_buffer.GetGenesBegin(idx);
    std::copy(genes.begin(), genes.end(), dest_begin);
    write_buffer.SetFitness(idx, fitness);
  }

  template <typename Iterator>
  void SetWriteIndividual(size_t idx, Iterator source_begin, double fitness) {
    auto dest_begin = write_buffer.GetGenesBegin(idx);
    std::copy(source_begin, source_begin + gene_size, dest_begin);
    write_buffer.SetFitness(idx, fitness);
  }

  [[nodiscard]]
  size_t GetGridSize() const {
    return read_buffer.GetGeneSize() * pop_size;
  }
  [[nodiscard]]
  size_t GetPopulationSize() const {
    return pop_size;
  }
  [[nodiscard]]
  size_t GetGeneSize() const {
    return read_buffer.GetGeneSize();
  }

  void FlushKBlocks(size_t start_idx, size_t k) {
    size_t buf_size = k * read_buffer.GetGeneSize();

    // Get memory blocks
    auto source_begin = write_buffer.GetGenesBegin(start_idx);
    auto source_end = source_begin + buf_size;
    auto dest_begin = read_buffer.GetGenesBegin(start_idx);

    // Copy block to read buffer
    std::copy(source_begin, source_end, dest_begin);
  }
};
