#include "afanasyev_a_integ_rect_method/seq/include/ops_seq.hpp"

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
  if (n <= 0) {
    return false;
  }

  const int k_dim = 3;

  const double h = 1.0 / static_cast<double>(n);

  std::vector<int> idx(k_dim, 0);
  std::vector<double> x(k_dim, 0.0);

  double sum = 0.0;

  bool done = false;
  while (!done) {
    for (int dim = 0; dim < k_dim; ++dim) {
      x[dim] = (static_cast<double>(idx[dim]) + 0.5) * h;
    }

    sum += ExampleIntegrand(x);

    for (int dim = 0; dim < k_dim; ++dim) {
      idx[dim]++;
      if (idx[dim] < n) {
        break;
      }
      idx[dim] = 0;
      if (dim == k_dim - 1) {
        done = true;
      }
    }
  }

  const double volume = std::pow(h, k_dim);
  GetOutput() = sum * volume;

  return true;
}

bool AfanasyevAIntegRectMethodSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace afanasyev_a_integ_rect_method
