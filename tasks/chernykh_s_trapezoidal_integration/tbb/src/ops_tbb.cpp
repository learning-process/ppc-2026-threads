#include "chernykh_s_trapezoidal_integration/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationTBB::ChernykhSTrapezoidalIntegrationTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationTBB::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationTBB::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationTBB::CalculatePointAndWeight(const IntegrationInType &input,
                                                                   const std::vector<std::size_t> &counters,
                                                                   std::vector<double> &point) {
  double weight = 1.0;
  for (std::size_t i = 0; i < input.limits.size(); ++i) {  // проходим по границам каждого из измерений
    const double h = (input.limits[i].second - input.limits[i].first) /
                     static_cast<double>(input.steps[i]);  // велечина шага в текущем измерении
    point[i] = input.limits[i].first +
               (static_cast<double>(counters[i]) * h);  // дискретная точка: начало_измерения + шаг*номер_шага
    if (std::cmp_equal(counters[i], 0) ||
        std::cmp_equal(counters[i], input.steps[i])) {  // если шаг граничный, то его вес уменьшается
      weight *= 0.5;
    }
  }
  return weight;
}

bool ChernykhSTrapezoidalIntegrationTBB::RunImpl() {
  const auto &input = this->GetInput();
  const std::size_t dims = input.limits.size();
  double total_sum = 0.0;
  int64_t total_points = 1;

  for (int steps_on_this_spase : input.steps) {
    total_points *= steps_on_this_spase + 1;
  }

  auto body = [&](const tbb::blocked_range<int64_t> &r, double local_sum) -> double {
    std::vector<double> local_point(dims);
    std::vector<size_t> local_counters(dims);

    for (int64_t j = r.begin(); j < r.end(); j++) {
      int64_t temp_j = j;
      for (size_t i = 0; i < input.steps.size(); i++) {
        local_counters[i] = temp_j % (input.steps[i] + 1);
        temp_j /= input.steps[i] + 1;
      }
      double weight = CalculatePointAndWeight(input, local_counters, local_point);
      local_sum += input.func(local_point) * weight;
    }

    return local_sum;
  };

  total_sum = tbb::parallel_reduce(tbb::blocked_range<int64_t>(0, total_points), 0.0, body, std::plus<>());

  double h_prod = 1.0;
  for (std::size_t i = 0; i < dims; ++i) {
    h_prod *= (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
  }

  GetOutput() = total_sum * h_prod;
  return true;
}

bool ChernykhSTrapezoidalIntegrationTBB::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
