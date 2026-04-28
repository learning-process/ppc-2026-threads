#include "shilin_n_monte_carlo_integration/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "shilin_n_monte_carlo_integration/common/include/common.hpp"

namespace shilin_n_monte_carlo_integration {

ShilinNMonteCarloIntegrationTBB::ShilinNMonteCarloIntegrationTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShilinNMonteCarloIntegrationTBB::ValidationImpl() {
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

bool ShilinNMonteCarloIntegrationTBB::PreProcessingImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  lower_bounds_ = lower;
  upper_bounds_ = upper;
  num_points_ = n;
  func_type_ = func_type;
  return true;
}

bool ShilinNMonteCarloIntegrationTBB::RunImpl() {
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

  double sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, num_points_), 0.0,
                                    [&](const tbb::blocked_range<int> &range, double local_sum) {
    std::vector<double> point(dimensions);
    for (int i = range.begin(); i < range.end(); ++i) {
      for (int di = 0; di < dimensions; ++di) {
        double val = 0.5 + (static_cast<double>(i + 1) * alpha[di]);
        double current = val - std::floor(val);
        point[di] = lower_bounds_[di] + ((upper_bounds_[di] - lower_bounds_[di]) * current);
      }
      local_sum += IntegrandFunction::Evaluate(func_type_, point);
    }
    return local_sum;
  }, std::plus<>());

  double volume = 1.0;
  for (int di = 0; di < dimensions; ++di) {
    volume *= (upper_bounds_[di] - lower_bounds_[di]);
  }

  GetOutput() = volume * sum / static_cast<double>(num_points_);
  return true;
}

bool ShilinNMonteCarloIntegrationTBB::PostProcessingImpl() {
  return true;
}

}  // namespace shilin_n_monte_carlo_integration
