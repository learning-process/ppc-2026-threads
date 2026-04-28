#include "redkina_a_integral_simpson/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"

namespace redkina_a_integral_simpson {

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

  const std::vector<double> a_local = a_;
  const std::vector<double> b_local = b_;
  const std::vector<int> n_local = n_;
  const auto func_local = func_;

  std::vector<double> h(dim);
  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_local[i] - a_local[i]) / static_cast<double>(n_local[i]);
    h_prod *= h[i];
  }

  std::vector<int> strides(dim);
  strides[dim - 1] = 1;
  for (int i = static_cast<int>(dim) - 2; i >= 0; --i) {
    strides[i] = strides[i + 1] * (n_local[i + 1] + 1);
  }
  int total_nodes = strides[0] * (n_local[0] + 1);

  double total_sum = 0.0;

#pragma omp parallel default(none) shared(total_nodes, h, strides, a_local, n_local, func_local, dim) \
    reduction(+ : total_sum)
  {
    std::vector<int> indices(dim);
    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (int idx = 0; idx < total_nodes; ++idx) {
      int remainder = idx;
      for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
        indices[dim_idx] = remainder / strides[dim_idx];
        remainder %= strides[dim_idx];
      }

      double w_prod = 1.0;
      for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
        int i_idx = indices[dim_idx];
        point[dim_idx] = a_local[dim_idx] + (i_idx * h[dim_idx]);

        int w = 0;
        if (i_idx == 0 || i_idx == n_local[dim_idx]) {
          w = 1;
        } else if (i_idx % 2 == 1) {
          w = 4;
        } else {
          w = 2;
        }
        w_prod *= static_cast<double>(w);
      }

      total_sum += w_prod * func_local(point);
    }
  }

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  result_ = (h_prod / denominator) * total_sum;
  return true;
}

bool RedkinaAIntegralSimpsonOMP::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson
