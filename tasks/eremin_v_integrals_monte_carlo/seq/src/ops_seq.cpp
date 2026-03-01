#include "eremin_v_integrals_monte_carlo/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <ranges>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"

namespace eremin_v_integrals_monte_carlo {

EreminVIntegralsMonteCarloSEQ::EreminVIntegralsMonteCarloSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool EreminVIntegralsMonteCarloSEQ::ValidationImpl() {
  const auto &input = GetInput();

  if (input.samples <= 0) {
    return false;
  }
  if (input.bounds.empty()) {
    return false;
  }
  if (!input.func) {
    return false;
  }

  return std::ranges::all_of(input.bounds, [](const auto &p) {
    const auto &[a, b] = p;
    return (a < b) && (std::abs(a) <= 1e9) && (std::abs(b) <= 1e9);
  });
  return true;
}

bool EreminVIntegralsMonteCarloSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool EreminVIntegralsMonteCarloSEQ::RunImpl() {
  const auto &input = GetInput();
  const auto &bounds = input.bounds;
  int samples = input.samples;
  const auto &func = input.func;

  const std::size_t dimension = bounds.size();

  double volume = 1.0;
  for (const auto &[a, b] : bounds) {
    volume *= (b - a);
  }

  std::mt19937 gen(std::random_device{}());

  std::vector<std::uniform_real_distribution<double>> distributions;
  distributions.reserve(dimension);

  for (const auto &[a, b] : bounds) {
    distributions.emplace_back(a, b);
  }

  double sum = 0.0;
  std::vector<double> point(dimension);

  for (int i = 0; i < samples; ++i) {
    for (std::size_t dim = 0; d < dimension; ++dim) {
      point[dim] = distributions[dim](gen);
    }

    sum += func(point);
  }

  GetOutput() = volume * (sum / static_cast<double>(samples));
  return true;
}

bool EreminVIntegralsMonteCarloSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace eremin_v_integrals_monte_carlo
