#include "afanasyev_a_integ_rect_method/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"

namespace afanasyev_a_integ_rect_method {
namespace {

double ExampleIntegrand(const std::array<double, 3> &x) {
  double s = 0.0;
  for (double xi : x) {
    s += xi * xi;
  }
  return std::exp(-s);
}

}  // namespace

AfanasyevAIntegRectMethodTBB::AfanasyevAIntegRectMethodTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool AfanasyevAIntegRectMethodTBB::ValidationImpl() {
  return (GetInput() > 0);
}

bool AfanasyevAIntegRectMethodTBB::PreProcessingImpl() {
  return true;
}

bool AfanasyevAIntegRectMethodTBB::RunImpl() {
  const int n = GetInput();
  if (n <= 0) {
    return false;
  }

  const int k_dim = 3;
  const double h = 1.0 / static_cast<double>(n);
  const int64_t total_points = static_cast<int64_t>(n) * n * n;

  const double sum = tbb::parallel_reduce(tbb::blocked_range<int64_t>(0, total_points), 0.0,
                                          [&](const tbb::blocked_range<int64_t> &range, double local_sum) {
    std::array<double, 3> x{};
    const int64_t plane = static_cast<int64_t>(n) * n;
    for (int64_t index = range.begin(); index != range.end(); ++index) {
      const int i = static_cast<int>(index / plane);
      const int j = static_cast<int>((index / n) % n);
      const int k = static_cast<int>(index % n);

      x[0] = (static_cast<double>(i) + 0.5) * h;
      x[1] = (static_cast<double>(j) + 0.5) * h;
      x[2] = (static_cast<double>(k) + 0.5) * h;

      local_sum += ExampleIntegrand(x);
    }
    return local_sum;
  }, std::plus<>());

  const double volume = std::pow(h, k_dim);
  GetOutput() = sum * volume;

  return true;
}

bool AfanasyevAIntegRectMethodTBB::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
