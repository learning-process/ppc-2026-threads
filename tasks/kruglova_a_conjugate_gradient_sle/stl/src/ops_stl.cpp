#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

double ParallelDotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
  return std::transform_reduce(std::execution::par, v1.begin(), v1.end(), v2.begin(), 0.0);
}

KruglovaAConjGradSleSTL::KruglovaAConjGradSleSTL(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleSTL::ValidationImpl() {
  const auto& in = GetInput();
  if (in.size <= 0) return false;
  if (in.A.size() != static_cast<size_t>(in.size) * in.size) return false;
  if (in.b.size() != static_cast<size_t>(in.size)) return false;
  return true;
}

bool KruglovaAConjGradSleSTL::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleSTL::RunImpl() {
  const auto& a = GetInput().A;
  const auto& b = GetInput().b;
  const int n = GetInput().size;
  auto& x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n);

  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  double rsold = ParallelDotProduct(r, r);
  const double tolerance = 1e-8;

  for (int iter = 0; iter < n * 2; ++iter) {

    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      double sum = 0.0;
      const size_t row_offset = static_cast<size_t>(i) * n;
      for (int j = 0; j < n; ++j) {
        sum += a[row_offset + j] * p[j];
      }
      ap[i] = sum;
    });

    double p_ap = ParallelDotProduct(p, ap);
    if (std::abs(p_ap) < 1e-15) break;

    const double alpha = rsold / p_ap;

    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    });

    const double rsnew = ParallelDotProduct(r, r);
    if (std::sqrt(rsnew) < tolerance) break;

    const double beta = rsnew / rsold;


    std::for_each(std::execution::par, indices.begin(), indices.end(), [&](int i) {
      p[i] = r[i] + beta * p[i];
    });

    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle