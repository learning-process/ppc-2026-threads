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
  // Защита от пустой функции (для анализатора)
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
  // Проверки для анализатора
  if (!func_) {
    return false;
  }
  size_t dim = a_.size();
  if (dim == 0) {
    return false;
  }

  // Шаги интегрирования по каждому измерению
  std::vector<double> h(dim);
  for (size_t i = 0; i < dim; ++i) {
    h[i] = (b_[i] - a_[i]) / static_cast<double>(n_[i]);
  }

  // Произведение шагов
  double h_prod = 1.0;
  for (size_t i = 0; i < dim; ++i) {
    h_prod *= h[i];
  }

  // Размеры сетки (количество точек в каждом измерении)
  std::vector<int> dim_sizes(dim);
  size_t total_points = 1;
  for (size_t i = 0; i < dim; ++i) {
    dim_sizes[i] = n_[i] + 1;
    total_points *= static_cast<size_t>(dim_sizes[i]);
  }

  // Вычисление strides (шагов для перехода между измерениями)
  std::vector<size_t> strides(dim);
  strides[dim - 1] = 1;
  for (size_t i = dim - 1; i > 0; --i) {
    strides[i - 1] = strides[i] * static_cast<size_t>(dim_sizes[i]);
  }

  // Параллельное суммирование с использованием std::async
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  size_t block_size = total_points / num_threads;
  size_t remainder_points = total_points % num_threads;

  std::vector<std::future<double>> futures;
  size_t current_start = 0;

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t current_end = current_start + block_size + (thread_idx < remainder_points ? 1 : 0);
    current_end = std::min(current_end, total_points);
    if (current_start >= current_end) {
      break;
    }

    // Запускаем асинхронную задачу. Захватываем необходимые данные:
    // - this: для доступа к a_, h_, n_, func_
    // - current_start, current_end: по значению
    // - strides: копируем, так как это локальный вектор
    futures.push_back(std::async(std::launch::async, [this, current_start, current_end, strides]() {
      double local_sum = 0.0;
      for (size_t idx = current_start; idx < current_end; ++idx) {
        local_sum += ComputeNodeContribution(idx, a_, h_, n_, strides, func_);
      }
      return local_sum;
    }));

    current_start = current_end;
  }

  double sum = 0.0;
  for (auto &fut : futures) {
    sum += fut.get();
  }

  // Знаменатель (3^dim)
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
