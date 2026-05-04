#include "telnov_a_integral_rectangle/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstdint>

#include "telnov_a_integral_rectangle/common/include/common.hpp"

namespace telnov_a_integral_rectangle {

TelnovAIntegralRectangleOMP::TelnovAIntegralRectangleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool TelnovAIntegralRectangleOMP::ValidationImpl() {
  return GetInput().first > 0 && GetInput().second > 0;
}

bool TelnovAIntegralRectangleOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TelnovAIntegralRectangleOMP::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  double result = 0.0;

#pragma omp parallel for default(none) reduction(+ : result) shared(total_points, n, d, a, h)
  for (int64_t idx = 0; idx < total_points; idx++) {
    int64_t tmp = idx;
    double f_value = 0.0;

    for (int dim = 0; dim < d; dim++) {
      const int coord_index = static_cast<int>(tmp % n);
      tmp /= n;

      const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
      f_value += x;
    }

    result += f_value;
  }

  result *= std::pow(h, d);

  GetOutput() = result;
  return true;
}

bool TelnovAIntegralRectangleOMP::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
