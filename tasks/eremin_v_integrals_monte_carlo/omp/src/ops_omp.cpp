#include "eremin_v_integrals_monte_carlo/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"

namespace eremin_v_integrals_monte_carlo {

EreminVIntegralsMonteCarloOMP::EreminVIntegralsMonteCarloOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool EreminVIntegralsMonteCarloOMP::ValidationImpl() {
  const auto &input = GetInput();

  if (input.samples <= 0) {
    return false;
  }
  if (input.bounds.empty()) {
    return false;
  }
  if (input.func == nullptr) {
    return false;
  }

  return std::ranges::all_of(input.bounds, [](const auto &p) {
    const auto &[a, b] = p;
    return (a < b) && (std::abs(a) <= 1e9) && (std::abs(b) <= 1e9);
  });
}

bool EreminVIntegralsMonteCarloOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool EreminVIntegralsMonteCarloOMP::RunImpl() {
  const auto &input = GetInput();
  const auto &bounds = input.bounds;
  int samples = input.samples;
  const auto &func = input.func;

  const std::size_t dimension = bounds.size();

  double volume = 1.0;
  for (const auto &[a, b] : bounds) {
    volume *= (b - a);
  }

  double sum = 0.0;
#pragma omp parallel reduction(+ : sum) default(none) shared(bounds, samples, func, dimension)
  {
    std::mt19937 local_gen(std::random_device{}() + omp_get_thread_num());

    std::vector<std::uniform_real_distribution<double>> local_distributions;
    local_distributions.reserve(dimension);

    for (const auto &[a, b] : bounds) {
      local_distributions.emplace_back(a, b);
    }

    std::vector<double> point(dimension);

#pragma omp for
    for (int i = 0; i < samples; ++i) {
      for (std::size_t dim = 0; dim < dimension; ++dim) {
        point[dim] = local_distributions[dim](local_gen);
      }

      sum += func(point);
    }
  }

  GetOutput() = volume * (sum / static_cast<double>(samples));
  return true;
}

bool EreminVIntegralsMonteCarloOMP::PostProcessingImpl() {
  return true;
}

}  // namespace eremin_v_integrals_monte_carlo
