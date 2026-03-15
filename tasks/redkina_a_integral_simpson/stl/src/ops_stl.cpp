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

// Вычисление вклада одного узла сетки
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

// Вычисление суммы на диапазоне индексов
double ComputeRange(size_t start, size_t end, const std::vector<double> &a, const std::vector<double> &h,
                    const std::vector<int> &n, const std::vector<size_t> &strides,
                    const std::function<double(const std::vector<double> &)> &func) {
  double local_sum = 0.0;
  for (size_t linear_idx = start; linear_idx < end; ++linear_idx) {
    local_sum += ComputeNodeContribution(linear_idx, a, h, n, strides, func);
  }
  return local_sum;
}

// Подготовка данных для интегрирования: шаги, размеры сетки, strides
bool PrepareData(const std::vector<double> &a, const std::vector<double> &b, const std::vector<int> &n,
                 std::vector<double> &h, std::vector<int> &dim_sizes, std::vector<size_t> &strides, double &h_prod,
                 size_t &total_points) {
  size_t dim = a.size();
  h.resize(dim);
  h_prod = 1.0;
  total_points = 1;
  dim_sizes.resize(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b[i] - a[i]) / static_cast<double>(n[i]);
    h_prod *= h[i];
    dim_sizes[i] = n[i] + 1;
    total_points *= static_cast<size_t>(dim_sizes[i]);
  }

  strides.resize(dim);
  strides[dim - 1] = 1;
  for (size_t i = dim - 1; i > 0; --i) {
    strides[i - 1] = strides[i] * static_cast<size_t>(dim_sizes[i]);
  }
  return true;
}

// Запуск параллельных задач и сбор результата
double ParallelSum(size_t total_points, const std::vector<double> &a, const std::vector<double> &h,
                   const std::vector<int> &n, const std::vector<size_t> &strides,
                   const std::function<double(const std::vector<double> &)> &func) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  size_t block_size = total_points / num_threads;
  size_t remainder_blocks = total_points % num_threads;

  std::vector<std::future<double>> futures;
  size_t start = 0;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t end = start + block_size + (thread_idx < remainder_blocks ? 1 : 0);
    end = std::min(end, total_points);
    if (start >= end) {
      break;
    }

    futures.push_back(std::async(std::launch::async, ComputeRange, start, end, std::cref(a), std::cref(h), std::cref(n),
                                 std::cref(strides), std::cref(func)));
    start = end;
    if (start >= total_points) {
      break;
    }
  }

  double sum = 0.0;
  for (auto &fut : futures) {
    sum += fut.get();
  }
  return sum;
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

  std::vector<double> h;
  std::vector<int> dim_sizes;
  std::vector<size_t> strides;
  double h_prod = 0.0;
  size_t total_points = 0;

  if (!PrepareData(a_, b_, n_, h, dim_sizes, strides, h_prod, total_points)) {
    return false;
  }

  double sum = ParallelSum(total_points, a_, h, n_, strides, func_);

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
