#include "shilin_n_monte_carlo_integration/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"

namespace shilin_n_monte_carlo_integration {

ShilinNMonteCarloIntegrationSTL::ShilinNMonteCarloIntegrationSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShilinNMonteCarloIntegrationSTL::ValidationImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  if (lower.size() != upper.size() || lower.empty()) {
    return false;
  }
  if (n <= 0) {
    return false;
  }
  for (size_t i = 0; i < lower.size(); ++i) {
    if (lower[i] >= upper[i]) {
      return false;
    }
  }
  if (func_type < FuncType::kConstant || func_type > FuncType::kSinProduct) {
    return false;
  }
  constexpr size_t kMaxDimensions = 10;
  return lower.size() <= kMaxDimensions;
}

bool ShilinNMonteCarloIntegrationSTL::PreProcessingImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  lower_bounds_ = lower;
  upper_bounds_ = upper;
  num_points_ = n;
  func_type_ = func_type;
  return true;
}

bool ShilinNMonteCarloIntegrationSTL::RunImpl() {
  auto dimensions = static_cast<int>(lower_bounds_.size());

  const std::vector<double> alpha = {
      0.41421356237309504,  // frac(sqrt(2))
      0.73205080756887729,  // frac(sqrt(3))
      0.23606797749978969,  // frac(sqrt(5))
      0.64575131106459059,  // frac(sqrt(7))
      0.31662479035539984,  // frac(sqrt(11))
      0.60555127546398929,  // frac(sqrt(13))
      0.12310562561766059,  // frac(sqrt(17))
      0.35889894354067355,  // frac(sqrt(19))
      0.79583152331271838,  // frac(sqrt(23))
      0.38516480713450403   // frac(sqrt(29))
  };

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<double> partial_sums(num_threads, 0.0);
  std::vector<std::thread> threads(num_threads);

  auto worker = [&](unsigned int tid) {
    std::vector<double> point(dimensions);
    for (int i = static_cast<int>(tid); i < num_points_; i += static_cast<int>(num_threads)) {
      for (int di = 0; di < dimensions; ++di) {
        double val = 0.5 + (static_cast<double>(i + 1) * alpha[di]);
        double current = val - std::floor(val);
        point[di] = lower_bounds_[di] + ((upper_bounds_[di] - lower_bounds_[di]) * current);
      }
      partial_sums[tid] += IntegrandFunction::Evaluate(func_type_, point);
    }
  };

  for (unsigned int ti = 0; ti < num_threads; ++ti) {
    threads[ti] = std::thread(worker, ti);
  }
  for (auto &th : threads) {
    th.join();
  }

  double sum = 0.0;
  for (unsigned int ti = 0; ti < num_threads; ++ti) {
    sum += partial_sums[ti];
  }

  double volume = 1.0;
  for (int di = 0; di < dimensions; ++di) {
    volume *= (upper_bounds_[di] - lower_bounds_[di]);
  }

  GetOutput() = volume * sum / static_cast<double>(num_points_);
  return true;
}

bool ShilinNMonteCarloIntegrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace shilin_n_monte_carlo_integration
