#include "telnov_a_integral_rectangle/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstdint>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "telnov_a_integral_rectangle/common/include/common.hpp"

namespace telnov_a_integral_rectangle {

TelnovAIntegralRectangleTBB::TelnovAIntegralRectangleTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool TelnovAIntegralRectangleTBB::ValidationImpl() {
  return GetInput().first > 0 && GetInput().second > 0;
}

bool TelnovAIntegralRectangleTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool TelnovAIntegralRectangleTBB::RunImpl() {
  const int n = GetInput().first;
  const int d = GetInput().second;

  const double a = 0.0;
  const double b = 1.0;
  const double h = (b - a) / static_cast<double>(n);

  const auto total_points = static_cast<int64_t>(std::pow(n, d));

  const double result =
      oneapi::tbb::parallel_reduce(oneapi::tbb::blocked_range<int64_t>(0, total_points), 0.0,
                                   [n, d, a, h](const oneapi::tbb::blocked_range<int64_t> &range, double local_sum) {
    for (int64_t idx = range.begin(); idx != range.end(); ++idx) {
      int64_t tmp = idx;
      double f_value = 0.0;

      for (int dim = 0; dim < d; ++dim) {
        const int coord_index = static_cast<int>(tmp % n);
        tmp /= n;

        const double x = a + ((static_cast<double>(coord_index) + 0.5) * h);
        f_value += x;
      }

      local_sum += f_value;
    }
    return local_sum;
  }, [](double lhs, double rhs) { return lhs + rhs; });

  GetOutput() = result * std::pow(h, d);
  return true;
}

bool TelnovAIntegralRectangleTBB::PostProcessingImpl() {
  return true;
}

}  // namespace telnov_a_integral_rectangle
