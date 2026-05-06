#include "redkina_a_integral_simpson/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"

namespace redkina_a_integral_simpson {

namespace {
double ComputeNodeContribution(size_t linear_idx, const std::vector<double> &a, const std::vector<double> &h,
                               const std::vector<int> &n, const std::vector<size_t> &strides,
                               const std::function<double(const std::vector<double> &)> &func) {
  size_t dim = a.size();
  std::vector<double> point(dim);
  size_t remainder = linear_idx;
  std::vector<int> indices(dim);

  for (size_t i = 0; i < dim; ++i) {
    indices[i] = static_cast<int>(remainder / strides[i]);
    remainder %= strides[i];
  }

  double w_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    int idx = indices[i];
    point[i] = a[i] + (static_cast<double>(idx) * h[i]);

    int w = 0;
    if (idx == 0 || idx == n[i]) {
      w = 1;
    } else if (idx % 2 == 1) {
      w = 4;
    } else {
      w = 2;
    }
    w_prod *= static_cast<double>(w);
  }
  return w_prod * func(point);
}
}  // namespace

RedkinaAIntegralSimpsonTBB::RedkinaAIntegralSimpsonTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonTBB::ValidationImpl() {
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

bool RedkinaAIntegralSimpsonTBB::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonTBB::RunImpl() {
  if (!func_) {
    return false;
  }
  size_t dim = a_.size();
  if (dim == 0) {
    return false;
  }

  std::vector<double> h(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
  }

  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h_prod *= h[i];
  }

  std::vector<int> dim_sizes(dim);
  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    dim_sizes[i] = n_[i] + 1;
    total_points *= static_cast<size_t>(dim_sizes[i]);
  }

  std::vector<size_t> strides(dim);
  strides[dim - 1] = 1;
  for (size_t i = dim - 1; i > 0; --i) {
    strides[i - 1] = strides[i] * static_cast<size_t>(dim_sizes[i]);
  }

  double sum = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, total_points), 0.0,
                                    [&](const tbb::blocked_range<size_t> &r, double local_sum) {
    for (size_t linear_idx = r.begin(); linear_idx != r.end(); ++linear_idx) {
      local_sum += ComputeNodeContribution(linear_idx, a_, h, n_, strides, func_);
    }
    return local_sum;
  }, [](double x, double y) -> double { return x + y; });

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  result_ = (h_prod / denominator) * sum;
  return true;
}

bool RedkinaAIntegralSimpsonTBB::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson
