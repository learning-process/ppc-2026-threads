#include "dergynov_s_integrals_multistep_rectangle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(borders, [](const auto &border) {
    return std::isfinite(border.first) && std::isfinite(border.second) && border.first < border.second;
  });
}

void ProcessRange(size_t start, size_t end, const std::function<double(const std::vector<double> &)> &func,
                  const std::vector<std::pair<double, double>> &borders, const std::vector<double> &h, int dim, int n,
                  std::atomic<bool> &error_flag, double &result) {
  double local_sum = 0.0;
  for (size_t linear_idx = start; linear_idx < end && !error_flag.load(); ++linear_idx) {
    size_t tmp = linear_idx;
    std::vector<double> point(dim);

    for (int dimension = dim - 1; dimension >= 0; --dimension) {
      int idx_val = static_cast<int>(tmp % static_cast<size_t>(n));
      tmp /= static_cast<size_t>(n);
      point[dimension] = borders[dimension].first + ((static_cast<double>(idx_val) + 0.5) * h[dimension]);
    }

    double f_val = func(point);
    if (!std::isfinite(f_val)) {
      error_flag.store(true);
      return;
    }
    local_sum += f_val;
  }
  result = local_sum;
}

}  // namespace

DergynovSIntegralsMultistepRectangleSTL::DergynovSIntegralsMultistepRectangleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleSTL::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();

  if (!func) {
    return false;
  }
  if (n <= 0) {
    return false;
  }
  if (borders.empty()) {
    return false;
  }

  return ValidateBorders(borders);
}

bool DergynovSIntegralsMultistepRectangleSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleSTL::RunImpl() {
  const auto &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;

  for (int idx = 0; idx < dim; ++idx) {
    const double left = borders[idx].first;
    const double right = borders[idx].second;
    h[idx] = (right - left) / static_cast<double>(n);
    cell_volume *= h[idx];
  }

  size_t total_points = 1;
  for (int idx = 0; idx < dim; ++idx) {
    total_points *= static_cast<size_t>(n);
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  std::vector<double> partial_sums(num_threads, 0.0);
  std::atomic<bool> error_flag{false};

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  size_t chunk_size = total_points / num_threads;
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start = thread_idx * chunk_size;
    size_t end = (thread_idx == num_threads - 1) ? total_points : start + chunk_size;

    threads.emplace_back(ProcessRange, start, end, func, borders, h, dim, n, std::ref(error_flag),
                         std::ref(partial_sums[thread_idx]));
  }

  for (auto &thread : threads) {
    thread.join();
  }

  if (error_flag.load()) {
    return false;
  }

  double total_sum = 0.0;
  for (double value : partial_sums) {
    total_sum += value;
  }

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleSTL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
