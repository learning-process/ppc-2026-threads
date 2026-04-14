#include "redkina_a_integral_simpson/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "redkina_a_integral_simpson/common/include/common.hpp"

namespace redkina_a_integral_simpson {

namespace {

double ComputeNodeContribution(size_t linear_idx, const std::vector<double> &a, const std::vector<double> &h,
                               const std::vector<int> &n, const std::vector<size_t> &strides,
                               const std::function<double(const std::vector<double> &)> &func) {
  if (!func) {
    return 0.0;
  }

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

double ParallelSum(const std::vector<double> &a, const std::vector<double> &h, const std::vector<int> &n,
                   const std::vector<size_t> &strides, const std::function<double(const std::vector<double> &)> &func,
                   size_t total_points) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  size_t block_size = total_points / num_threads;
  size_t remainder = total_points % num_threads;

  std::vector<std::future<double>> futures;
  size_t start = 0;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t end = std::min(start + block_size + (thread_idx < remainder ? 1 : 0), total_points);
    if (start >= end) {
      break;
    }

    futures.push_back(std::async(std::launch::async, [=, &a, &h, &n, &strides, &func]() {
      double local_sum = 0.0;
      for (size_t idx = start; idx < end; ++idx) {
        local_sum += ComputeNodeContribution(idx, a, h, n, strides, func);
      }
      return local_sum;
    }));

    start = end;
    if (start >= total_points) {
      break;
    }
  }

  double total = 0.0;
  for (auto &f : futures) {
    total += f.get();
  }
  return total;
}

}  // namespace

RedkinaAIntegralSimpsonSTL::RedkinaAIntegralSimpsonSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RedkinaAIntegralSimpsonSTL::ValidationImpl() {
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

bool RedkinaAIntegralSimpsonSTL::PreProcessingImpl() {
  const auto &in = GetInput();
  func_ = in.func;
  a_ = in.a;
  b_ = in.b;
  n_ = in.n;
  result_ = 0.0;
  return true;
}

bool RedkinaAIntegralSimpsonSTL::RunImpl() {
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

  double sum = ParallelSum(a_, h, n_, strides, func_, total_points);

  double denominator = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    denominator *= 3.0;
  }

  result_ = (h_prod / denominator) * sum;
  return true;
}

bool RedkinaAIntegralSimpsonSTL::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace redkina_a_integral_simpson
