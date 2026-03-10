#include "redkina_a_integral_simpson_seq/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"

namespace redkina_a_integral_simpson_seq {

namespace {

inline int SimpsonCoeff(int idx, int n) {
  if (idx == 0 || idx == n) {
    return 1;
  }
  return (idx % 2 == 1) ? 4 : 2;
}

bool AdvanceIndicesFromLevel(std::vector<int> &indices, const std::vector<int> &n, int start_level) {
  int dim = static_cast<int>(indices.size());
  int d = dim - 1;
  while (d >= start_level && indices[d] == n[d]) {
    indices[d] = 0;
    --d;
  }
  if (d < start_level) {
    return false;
  }
  ++indices[d];
  return true;
}

}  // namespace

RedkinaAIntegralSimpsonOMP::RedkinaAIntegralSimpsonOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonOMP::ValidationImpl() {
  const auto &in = GetInput();
  size_t dim = in.a.size();

  if (dim == 0 || in.b.size() != dim || in.n.size() != dim) {
    return false;
  }

  for (size_t i = 0; i < dim; ++i) {
    if (in.a[i] >= in.b[i]) {
      return false;
    }
    if (in.n[i] <= 0 || in.n[i] % 2 != 0) {
      return false;
    }
  }

  return static_cast<bool>(in.func);
}

bool RedkinaAIntegralSimpsonOMP::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonOMP::RunImpl() {
  size_t dim = a_.size();

  std::vector<double> h(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
  }

  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h_prod *= h[i];
  }

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  // Локальные переменные для использования в OpenMP
  const std::vector<double> &a_ref = a_;
  const std::vector<double> &h_ref = h;
  const std::vector<int> &n_ref = n_;
  const auto &func_ref = func_;
  const size_t dim_local = dim;

  double total_sum = 0.0;

#pragma omp parallel for reduction(+ : total_sum) default(none) shared(a_ref, h_ref, n_ref, func_ref, dim_local)
  for (int i0 = 0; i0 <= static_cast<int>(n_ref[0]); ++i0) {
    double coeff0 = SimpsonCoeff(i0, static_cast<int>(n_ref[0]));
    double local_sum = 0.0;

    std::vector<double> point(dim_local);
    std::vector<int> indices(dim_local);
    indices[0] = i0;
    for (size_t d = 1; d < dim_local; ++d) {
      indices[d] = 0;
    }

    do {
      // Координата первого измерения
      point[0] = a_ref[0] + static_cast<double>(i0) * h_ref[0];

      double w_prod = 1.0;
      // Остальные измерения
      for (size_t d = 1; d < dim_local; ++d) {
        int idx = indices[d];
        point[d] = a_ref[d] + static_cast<double>(idx) * h_ref[d];
        int w = SimpsonCoeff(idx, static_cast<int>(n_ref[d]));
        w_prod *= static_cast<double>(w);
      }

      local_sum += coeff0 * w_prod * func_ref(point);
    } while (AdvanceIndicesFromLevel(indices, n_ref, 1));

    total_sum += local_sum;
  }

  result_ = (h_prod / denominator) * total_sum;
  return true;
}

bool RedkinaAIntegralSimpsonOMP::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson_seq
