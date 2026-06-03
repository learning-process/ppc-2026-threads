#include "eremin_v_integrals_monte_carlo/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <thread>
#include <vector>

#include "eremin_v_integrals_monte_carlo/common/include/common.hpp"
#include "util/include/util.hpp"

namespace eremin_v_integrals_monte_carlo {

EreminVIntegralsMonteCarloSTL::EreminVIntegralsMonteCarloSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool EreminVIntegralsMonteCarloSTL::ValidationImpl() {
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

bool EreminVIntegralsMonteCarloSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool EreminVIntegralsMonteCarloSTL::RunImpl() {
  const auto &input = GetInput();
  const auto &bounds = input.bounds;
  const int samples = input.samples;
  const auto &func = input.func;

  const std::size_t dimension = bounds.size();

  double volume = 1.0;
  for (const auto &[a, b] : bounds) {
    volume *= (b - a);
  }

  const int num_threads = std::max(1, ppc::util::GetNumThreads());
  std::vector<double> partial_sums(static_cast<std::size_t>(num_threads), 0.0);
  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(num_threads));

  const unsigned base_seed = std::random_device{}();

  for (int tid = 0; tid < num_threads; ++tid) {
    const int begin = (samples * tid) / num_threads;
    const int end = (samples * (tid + 1)) / num_threads;
    threads.emplace_back([&, tid, begin, end]() {
      std::minstd_rand local_gen(base_seed + static_cast<unsigned>(tid));
      std::vector<std::uniform_real_distribution<double>> local_distributions;
      local_distributions.reserve(dimension);
      for (const auto &[a, b] : bounds) {
        local_distributions.emplace_back(a, b);
      }

      std::vector<double> point(dimension);
      double local_sum = 0.0;
      for (int i = begin; i < end; ++i) {
        for (std::size_t dim = 0; dim < dimension; ++dim) {
          point[dim] = local_distributions[dim](local_gen);
        }
        local_sum += func(point);
      }
      partial_sums[static_cast<std::size_t>(tid)] = local_sum;
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  double sum = 0.0;
  for (double s : partial_sums) {
    sum += s;
  }

  GetOutput() = volume * (sum / static_cast<double>(samples));
  return true;
}

bool EreminVIntegralsMonteCarloSTL::PostProcessingImpl() {
  return true;
}

}  // namespace eremin_v_integrals_monte_carlo
