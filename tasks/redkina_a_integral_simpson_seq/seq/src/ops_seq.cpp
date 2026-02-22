#include "redkina_a_integral_simpson_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "redkina_a_integral_simpson_seq/common/include/common.hpp"

namespace redkina_a_integral_simpson_seq {

namespace {

void EvaluatePoint(const std::vector<double> &a, const std::vector<double> &h, const std::vector<int> &n,
                   const std::vector<int> &indices, const std::function<double(const std::vector<double> &)> &func,
                   std::vector<double> &point, double &sum) {
  size_t dim = a.size();
  double w_prod = 1.0;
  for (size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
    int idx = indices[dim_idx];
    point[dim_idx] = a[dim_idx] + (static_cast<double>(idx) * h[dim_idx]);

    int w = 0;
    if (idx == 0 || idx == n[dim_idx]) {
      w = 1;
    } else if (idx % 2 == 1) {
      w = 4;
    } else {
      w = 2;
    }
    w_prod *= static_cast<double>(w);
  }
  sum += w_prod * func(point);
}

bool AdvanceIndices(std::vector<int> &indices, const std::vector<int> &n) {
  int dim = static_cast<int>(indices.size());
  int d = dim - 1;
  while (d >= 0 && indices[d] == n[d]) {
    indices[d] = 0;
    --d;
  }
  if (d < 0) {
    return false;
  }
  ++indices[d];
  return true;
}

}  // namespace

RedkinaAIntegralSimpsonSEQ::RedkinaAIntegralSimpsonSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonSEQ::ValidationImpl() {
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

bool RedkinaAIntegralSimpsonSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonSEQ::RunImpl() {
  size_t dim = a_.size();

  std::vector<double> h(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
  }

  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h_prod *= h[i];
  }

  std::vector<double> point(dim);
  double sum = 0.0;

  std::vector<int> indices(dim, 0);

  bool has_next = true;
  while (has_next) {
    EvaluatePoint(a_, h, n_, indices, func_, point, sum);
    has_next = AdvanceIndices(indices, n_);
  }

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  result_ = (h_prod / denominator) * sum;
  return true;
}

bool RedkinaAIntegralSimpsonSEQ::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson_seq
