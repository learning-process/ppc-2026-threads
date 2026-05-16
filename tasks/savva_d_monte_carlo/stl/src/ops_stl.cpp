#include "savva_d_monte_carlo/stl/include/ops_stl.hpp"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"

namespace savva_d_monte_carlo {

SavvaDMonteCarloSTL::SavvaDMonteCarloSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SavvaDMonteCarloSTL::ValidationImpl() {
  const auto &input = GetInput();

  if (input.count_points == 0) {
    return false;
  }
  if (!input.f) {
    return false;
  }
  if (input.Dimension() == 0) {
    return false;
  }

  for (size_t i = 0; i < input.Dimension(); ++i) {
    if (input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
  }

  return true;
}

bool SavvaDMonteCarloSTL::PreProcessingImpl() {
  return true;
}

bool SavvaDMonteCarloSTL::RunImpl() {
  const auto &input = GetInput();
  auto &result = GetOutput();

  const size_t dim = input.Dimension();
  const double vol = input.Volume();
  const auto n = static_cast<int64_t>(input.count_points);
  const auto &func = input.f;

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  if (static_cast<unsigned int>(n) < num_threads) {
    num_threads = static_cast<unsigned int>(n);
  }

  std::vector<double> partial_sums(num_threads, 0.0);
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  const int64_t points_per_thread = n / num_threads;
  const int64_t tail = n % num_threads;

  auto worker = [&](int thread_id, int64_t pts) {
    std::minstd_rand gen(1337 + thread_id);

    std::vector<std::uniform_real_distribution<double>> dists;
    dists.reserve(dim);
    for (size_t i = 0; i < dim; ++i) {
      dists.emplace_back(input.lower_bounds[i], input.upper_bounds[i]);
    }

    std::vector<double> point(dim);

    double local_sum = 0.0;

    for (int64_t i = 0; i < pts; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        point[j] = dists[j](gen);
      }
      local_sum += func(point);
    }

    partial_sums[thread_id] = local_sum;
  };

  for (unsigned int i = 0; i < num_threads; ++i) {
    int64_t pts = points_per_thread + (i < static_cast<unsigned int>(tail) ? 1 : 0);
    threads.emplace_back(worker, i, pts);
  }

  for (auto &t : threads) {
    t.join();
  }

  double total_sum = 0.0;
  for (double s : partial_sums) {
    total_sum += s;
  }

  result = vol * total_sum / static_cast<double>(n);
  return true;
}

bool SavvaDMonteCarloSTL::PostProcessingImpl() {
  return true;
}

}  // namespace savva_d_monte_carlo
