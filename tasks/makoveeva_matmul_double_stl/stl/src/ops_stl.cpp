#include "makoveeva_matmul_double_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "makoveeva_matmul_double_stl/common/include/common.hpp"

namespace makoveeva_matmul_double_stl {
namespace {

void MultiplyRowRange(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n,
                      size_t start_row, size_t end_row) {
  for (size_t i = start_row; i < end_row; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] = sum;
    }
  }
}

}  // namespace

MatmulDoubleSTLTask::MatmulDoubleSTLTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleSTLTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleSTLTask::PreProcessingImpl() {
  const auto &input = GetInput();
  n_ = std::get<0>(input);
  A_ = std::get<1>(input);
  B_ = std::get<2>(input);
  C_.assign(n_ * n_, 0.0);

  return true;
}

bool MatmulDoubleSTLTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  const size_t num_threads = std::thread::hardware_concurrency();
  const size_t rows_per_thread = (n + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;

  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start_row = thread_idx * rows_per_thread;
    size_t end_row = (std::min)(start_row + rows_per_thread, n);

    if (start_row >= n) {
      break;
    }

    threads.emplace_back([&a, &b, &c, n, start_row, end_row]() { MultiplyRowRange(a, b, c, n, start_row, end_row); });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  GetOutput() = C_;
  return true;
}

bool MatmulDoubleSTLTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_stl
