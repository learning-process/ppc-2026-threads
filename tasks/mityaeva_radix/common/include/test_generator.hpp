#pragma once

#include <cstddef>
#include <random>
#include <vector>

namespace mityaeva_radix {

inline std::vector<double> GenerateTest(std::size_t length, std::size_t seed) {
  std::vector<double> input;
  input.reserve(length);
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (size_t i = 0; i < length; i++) {
    input.push_back(dist(rng) - 0.5);
  }

  return input;
};

}  // namespace mityaeva_radix
