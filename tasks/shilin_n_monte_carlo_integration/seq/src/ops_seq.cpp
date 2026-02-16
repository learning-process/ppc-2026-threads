#include "shilin_n_monte_carlo_integration/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "shilin_n_monte_carlo_integration/common/include/common.hpp"

namespace shilin_n_monte_carlo_integration {

ShilinNMonteCarloIntegrationSEQ::ShilinNMonteCarloIntegrationSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShilinNMonteCarloIntegrationSEQ::ValidationImpl() {
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
  return true;
}

bool ShilinNMonteCarloIntegrationSEQ::PreProcessingImpl() {
  const auto &[lower, upper, n, func_type] = GetInput();
  lower_bounds_ = lower;
  upper_bounds_ = upper;
  num_points_ = n;
  func_type_ = func_type;
  return true;
}

bool ShilinNMonteCarloIntegrationSEQ::RunImpl() {
  auto dimensions = static_cast<int>(lower_bounds_.size());

  // Quasi-random Kronecker sequence: alpha_d = fractional part of sqrt(prime_d)
  const std::vector<double> kAlpha = {
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

  auto alpha_size = static_cast<int>(kAlpha.size());
  std::vector<double> current(dimensions, 0.5);
  std::vector<double> point(dimensions);
  double sum = 0.0;

  for (int i = 0; i < num_points_; ++i) {
    for (int d = 0; d < dimensions; ++d) {
      current[d] += kAlpha[d % alpha_size];
      if (current[d] >= 1.0) {
        current[d] -= 1.0;
      }
      point[d] = lower_bounds_[d] + (upper_bounds_[d] - lower_bounds_[d]) * current[d];
    }
    sum += IntegrandFunction::Evaluate(func_type_, point);
  }

  double volume = 1.0;
  for (int d = 0; d < dimensions; ++d) {
    volume *= (upper_bounds_[d] - lower_bounds_[d]);
  }

  GetOutput() = volume * sum / static_cast<double>(num_points_);
  return true;
}

bool ShilinNMonteCarloIntegrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shilin_n_monte_carlo_integration
