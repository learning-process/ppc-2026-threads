#include "telnov_a_integral_rectangle/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "telnov_a_integral_rectangle/common/include/common.hpp"
#include "util/include/util.hpp"

namespace telnov_a_integral_rectangle {

TelnovAIntegralRectangleSEQ::TelnovAIntegralRectangleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool TelnovAIntegralRectangleSEQ::ValidationImpl() {
  return GetInput().first > 0 && GetInput().second > 0;
}

bool TelnovAIntegralRectangleSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TelnovAIntegralRectangleSEQ::RunImpl() {
  const int N = GetInput().first;
  const int D = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / N;

  const long long total_points = static_cast<long long>(std::pow(N, D));

  double result = 0.0;

  for (long long idx = 0; idx < total_points; idx++) {
    long long tmp = idx;
    double f_value = 0.0;

    for (int dim = 0; dim < D; dim++) {
      int coord_index = tmp % N;
      tmp /= N;

      double x = a + (coord_index + 0.5) * h;
      f_value += x;
    }

    result += f_value;
  }

  result *= std::pow(h, D);

  GetOutput() = result;
  return true;
}

bool TelnovAIntegralRectangleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
