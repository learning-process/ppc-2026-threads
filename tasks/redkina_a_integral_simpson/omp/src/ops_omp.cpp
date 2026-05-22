#include "redkina_a_integral_simpson/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"

namespace redkina_a_integral_simpson {
namespace {

std::vector<std::vector<double>> PrecomputeWeights(const std::vector<int> &n) {
  const size_t dim = n.size();
  std::vector<std::vector<double>> weights(dim);
  for (size_t i = 0; i < dim; ++i) {
    const int ni = n[i];
    weights[i].resize(ni + 1);
    for (int idx = 0; idx <= ni; ++idx) {
      if (idx == 0 || idx == ni) {
        weights[i][idx] = 1.0;
      } else if (idx % 2 == 1) {
        weights[i][idx] = 4.0;
      } else {
        weights[i][idx] = 2.0;
      }
    }
  }
  return weights;
}

std::vector<size_t> ComputeStrides(const std::vector<int> &n) {
  const size_t dim = n.size();
  std::vector<size_t> strides(dim);
  if (dim == 0) {
    return strides;
  }
  strides[dim - 1] = 1;
  for (size_t i = dim - 1; i > 0; --i) {
    strides[i - 1] = strides[i] * static_cast<size_t>(n[i] + 1);
  }
  return strides;
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
  if (!func_) {
    return false;
  }
  const size_t dim = a_.size();
  if (dim == 0) {
    return false;
  }

  std::vector<double> h(dim);
  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
    h_prod *= h[i];
  }

  const auto weights = PrecomputeWeights(n_);
  const auto strides = ComputeStrides(n_);
  if (strides.empty()) {
    return false;
  }

  const int total_nodes = static_cast<int>(strides[0] * static_cast<size_t>(n_[0] + 1));

  double total_sum = 0.0;

  const std::vector<double> &a_local = a_;
  const std::vector<double> &h_local = h;
  const auto &weights_local = weights;
  const auto &strides_local = strides;
  const auto &func_local = func_;

#pragma omp parallel default(none) \
    shared(total_nodes, a_local, h_local, weights_local, strides_local, func_local, dim) reduction(+ : total_sum)
  {
    std::vector<int> indices(dim);
    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (int idx = 0; idx < total_nodes; ++idx) {
      auto remainder = static_cast<size_t>(idx);
      for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
        indices[dim_idx] = static_cast<int>(remainder / strides_local[dim_idx]);
        remainder %= strides_local[dim_idx];
      }

      double w_prod = 1.0;
      for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
        const int i_idx = indices[dim_idx];
        point[dim_idx] = a_local[dim_idx] + (static_cast<double>(i_idx) * h_local[dim_idx]);
        w_prod *= weights_local[dim_idx][i_idx];
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
