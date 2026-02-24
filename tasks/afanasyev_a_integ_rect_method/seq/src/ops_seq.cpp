#include "afanasyev_a_integ_rect_method/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>
#include <functional>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"

namespace afanasyev_a_integ_rect_method {

static double ExampleIntegrand(const std::vector<double> &x) {
  double s = 0.0;
  for (double xi : x) s += xi * xi;
  return std::exp(-s);
}

AfanasyevAIntegRectMethodSEQ::AfanasyevAIntegRectMethodSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool AfanasyevAIntegRectMethodSEQ::ValidationImpl() {
  return (GetInput() > 0);
}

bool AfanasyevAIntegRectMethodSEQ::PreProcessingImpl() {
  return true;
}

bool AfanasyevAIntegRectMethodSEQ::RunImpl() {
  const int n = GetInput();
  if (n <= 0) return false;

  const int kDim = 3; 

  const double h = 1.0 / static_cast<double>(n);

  std::vector<int> idx(kDim, 0);
  std::vector<double> x(kDim, 0.0);

  double sum = 0.0;

  bool done = false;
  while (!done) {
    for (int d = 0; d < kDim; ++d) {
      x[d] = (static_cast<double>(idx[d]) + 0.5) * h;
    }

    sum += ExampleIntegrand(x);

    for (int d = 0; d < kDim; ++d) {
      idx[d]++;
      if (idx[d] < n) break;
      idx[d] = 0;
      if (d == kDim - 1) done = true;
    }
  }

  const double volume = std::pow(h, kDim);
  GetOutput() = sum * volume;

  return true;
}

bool AfanasyevAIntegRectMethodSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
