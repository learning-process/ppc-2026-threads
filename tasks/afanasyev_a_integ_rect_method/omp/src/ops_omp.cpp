#include "afanasyev_a_integ_rect_method/omp/include/ops_omp.hpp"

#include <cmath>
#include <vector>

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"

namespace afanasyev_a_integ_rect_method {
namespace {

double ExampleIntegrand(const std::vector<double> &x) {
  double s = 0.0;
  for (double xi : x) {
    s += xi * xi;
  }
  return std::exp(-s);
}

}  // namespace

AfanasyevAIntegRectMethodOMP::AfanasyevAIntegRectMethodOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool AfanasyevAIntegRectMethodOMP::ValidationImpl() {
  return (GetInput() > 0);
}

bool AfanasyevAIntegRectMethodOMP::PreProcessingImpl() {
  return true;
}

bool AfanasyevAIntegRectMethodOMP::RunImpl() {
  const int n = GetInput();
  if (n <= 0) {
    return false;
  }

  const int k_dim = 3;

  const double h = 1.0 / static_cast<double>(n);

  double sum = 0.0;

#pragma omp parallel for collapse(3) reduction(+ : sum) default(none) shared(n, h, k_dim)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        std::vector<double> x = {(static_cast<double>(i) + 0.5) * h, (static_cast<double>(j) + 0.5) * h,
                                 (static_cast<double>(k) + 0.5) * h};
        sum += ExampleIntegrand(x);
      }
    }
  }

  const double volume = std::pow(h, k_dim);
  GetOutput() = sum * volume;

  return true;
}

bool AfanasyevAIntegRectMethodOMP::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
