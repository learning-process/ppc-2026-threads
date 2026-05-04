#include "eremin_v_integrals_monte_carlo/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/parallel_reduce.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <random>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"

namespace eremin_v_integrals_monte_carlo {

EreminVIntegralsMonteCarloTBB::EreminVIntegralsMonteCarloTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool EreminVIntegralsMonteCarloTBB::ValidationImpl() {
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

bool EreminVIntegralsMonteCarloTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool EreminVIntegralsMonteCarloTBB::RunImpl() {
  const auto &input = GetInput();
  const auto &bounds = input.bounds;
  const int samples = input.samples;
  const auto &func = input.func;

  const std::size_t dimension = bounds.size();

  double volume = 1.0;
  for (const auto &[a, b] : bounds) {
    volume *= (b - a);
  }

  const double sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, samples), 0.0,
                                          [&](const tbb::blocked_range<int> &range, double local_sum) {
    // Thread-local RNG and distributions for this range chunk.
    std::mt19937 local_gen(std::random_device{}() + static_cast<unsigned>(range.begin()));

    std::vector<std::uniform_real_distribution<double>> local_distributions;
    local_distributions.reserve(dimension);
    for (const auto &[a, b] : bounds) {
      local_distributions.emplace_back(a, b);
    }

    std::vector<double> point(dimension);
    for (int i = range.begin(); i < range.end(); ++i) {
      for (std::size_t dim = 0; dim < dimension; ++dim) {
        point[dim] = local_distributions[dim](local_gen);
      }
      local_sum += func(point);
    }
    return local_sum;
  }, std::plus<double>{});

  GetOutput() = volume * (sum / static_cast<double>(samples));
  return true;
}

bool EreminVIntegralsMonteCarloTBB::PostProcessingImpl() {
  return true;
}

}  // namespace eremin_v_integrals_monte_carlo
