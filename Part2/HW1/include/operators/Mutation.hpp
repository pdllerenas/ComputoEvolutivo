#pragma once

#include <random>
#include <vector>
#include <cassert>

template <typename T>
void UniformMutation(std::vector<T> &child, double mutation_rate, T lower, T upper, std::mt19937 &rng) {
	assert(mutation_rate <= 1.0 && mutation_rate >= 0.0);
	assert(lower <= upper);
	
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	std::uniform_int_distribution<T> val_dis(lower, upper);
	for (size_t i = 0; i < child.size(); ++i) {
		if (dis(rng) <= mutation_rate) {
			child[i] = val_dis(rng);
		}
	}
}
