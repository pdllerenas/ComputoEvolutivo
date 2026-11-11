#pragma once

#include <iterator>
#include <random>
#include <algorithm>

template <typename Iterator,
          typename T = typename std::iterator_traits<Iterator>::value_type>
std::vector<T> BoundedDiscreteCrossover(Iterator it_mother, Iterator it_father,
                                        size_t dimension, std::mt19937 &rng) {

  std::vector<T> child(dimension);

  for (size_t i = 0; i < dimension; ++i) {
    auto [min_val, max_val] = std::minmax(*it_mother, *it_father);

    std::uniform_int_distribution<T> dis(min_val, max_val);
    child[i] = dis(rng);

    it_mother++;
    it_father++;
  }
  return child;
}
