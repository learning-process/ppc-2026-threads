#include "chernykh_s_trapezoidal_integration/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationSEQ::ChernykhSTrapezoidalIntegrationSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationSEQ::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationSEQ::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationSEQ::CalculatePointAndWeight(const IntegrationInType &input,
                                                                   const std::vector<std::size_t> &counters,
                                                                   std::vector<double> &point) {
  double weight = 1.0;
  for (std::size_t i = 0; i < input.limits.size(); ++i) {
    const double h = (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
    point[i] = input.limits[i].first + (static_cast<double>(counters[i]) * h);
    if (std::cmp_equal(counters[i], 0) || std::cmp_equal(counters[i], input.steps[i])) {
      weight *= 0.5;
    }
  }
  return weight;
}

bool ChernykhSTrapezoidalIntegrationSEQ::RunImpl() {
  const auto &input = this->GetInput();
  const std::size_t dims = input.limits.size();
  std::vector<std::size_t> counters(dims, 0);
  std::vector<double> current_point(dims);
  double total_sum = 0.0;
  bool done = false;

  while (!done) {
    double weight = CalculatePointAndWeight(input, counters, current_point);
    total_sum += input.func(current_point) * weight;

    for (std::size_t i = 0; i < dims; ++i) {
      if (std::cmp_less(++counters[i], input.steps[i] + 1)) {
        break;
      }
      if (std::cmp_equal(i, dims - 1)) {
        done = true;
      } else {
        counters[i] = 0;
      }
    }
  }

  double h_prod = 1.0;
  for (std::size_t i = 0; i < dims; ++i) {
    h_prod *= (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
  }

  GetOutput() = total_sum * h_prod;
  return true;
}

bool ChernykhSTrapezoidalIntegrationSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
