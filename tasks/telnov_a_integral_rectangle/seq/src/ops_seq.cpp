#include "telnov_a_integral_rectangle/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstdint>

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
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / n;

  const int64_t totalPoints = static_cast<int64_t>(std::pow(n, d));

  double result = 0.0;

  for (int64_t idx = 0; idx < totalPoints; idx++) {
    int64_t tmp = idx;
    double fValue = 0.0;

    for (int dim = 0; dim < d; dim++) {
      int coordIndex = static_cast<int>(tmp % n);
      tmp /= n;

      double x = a + ((coordIndex + 0.5) * h);
      fValue += x;
    }

    result += fValue;
  }

  result *= std::pow(h, d);

  GetOutput() = result;
  return true;
}

bool TelnovAIntegralRectangleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
