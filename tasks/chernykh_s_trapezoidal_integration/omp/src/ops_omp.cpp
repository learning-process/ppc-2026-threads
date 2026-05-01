#include "chernykh_s_trapezoidal_integration/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationOMP::ChernykhSTrapezoidalIntegrationOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationOMP::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationOMP::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationOMP::CalculatePointAndWeight(const IntegrationInType &input,
                                                                   const std::vector<std::size_t> &counters,
                                                                   std::vector<double> &point) {
  double weight = 1.0;
  for (std::size_t i = 0; i < input.limits.size(); ++i) {
    const double h = (input.limits[i].second - input.limits[i].first) /
                     static_cast<double>(input.steps[i]);  // шаг сетки h по i-ому измерению
    point[i] =
        input.limits[i].first + (static_cast<double>(counters[i]) * h);  // координата текущей точки в i измерении
    if (std::cmp_equal(counters[i], 0) ||
        std::cmp_equal(counters[i], input.steps[i])) {  // если это граничная точка, уменьшаем вес на половину
      weight *= 0.5;
    }
  }
  return weight;
}

bool ChernykhSTrapezoidalIntegrationOMP::RunImpl() {
  const auto &input = this->GetInput();
  const std::size_t dims = input.limits.size();
  int64_t total_points = 1;
  for (int setka : input.steps) {
    total_points *= (static_cast<int64_t>(setka) + 1);  // растет очень быстро
  }

  double total_sum = 0.0;
#pragma omp parallel default(none) shared(input, dims, total_points) reduction(+ : total_sum)
  {
    std::vector<std::size_t> local_counters(dims);  // создаем локальный вектор итераций
    std::vector<double> local_point(dims);          // значения в точках

#pragma omp for schedule(static)
    for (int64_t j = 0; j < total_points; j++) {
      int64_t temp_j = j;
      for (std::size_t i = 0; i < dims; i++) {
        int64_t point_in_dims = static_cast<int64_t>(input.steps[i]) + 1;  // количество точек в текущем измерении
        local_counters[i] = static_cast<std::size_t>(temp_j % point_in_dims);
        temp_j /= point_in_dims;
      }
      double weight = CalculatePointAndWeight(input, local_counters, local_point);
      total_sum += input.func(local_point) * weight;
    }
  }

  double h_prod = 1.0;
  for (std::size_t i = 0; i < dims; ++i) {
    h_prod *= (input.limits[i].second - input.limits[i].first) / static_cast<double>(input.steps[i]);
  }

  GetOutput() = total_sum * h_prod;
  return true;
}

bool ChernykhSTrapezoidalIntegrationOMP::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
