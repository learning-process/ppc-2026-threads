#include "makoveeva_matmul_double/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include "makoveeva_matmul_double/stl/include/common.hpp"

namespace makoveeva_matmul_double_stl {

namespace {

[[nodiscard]] size_t SelectBlockSize(size_t n) {
  if (n <= 64) {
    return n;
  }
  if (n <= 256) {
    return 64;
  }
  if (n <= 1024) {
    return 128;
  }
  return 256;
}
void SimpleMultiplyThread(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n,
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

  const size_t block_size = SelectBlockSize(n);

  if (n % block_size != 0) {
    return RunSimpleMultiply();
  }

  const size_t grid_size = n / block_size;

  const size_t num_threads = std::thread::hardware_concurrency();

  std::mutex write_mutex;

  const size_t total_iterations = grid_size * grid_size * grid_size;

  if (total_iterations >= num_threads) {
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    const size_t iterations_per_thread = total_iterations / num_threads;

    // Создаём потоки
    for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      const size_t start_step = thread_idx * iterations_per_thread;
      const size_t end_step = (thread_idx == num_threads - 1) ? total_iterations : start_step + iterations_per_thread;

      threads.emplace_back(&MatmulDoubleSTLTask::Worker, this, start_step, end_step, grid_size, block_size,
                           std::ref(write_mutex));
    }

    for (auto &thread : threads) {
      thread.join();
    }
  } else {
    Worker(0, total_iterations, grid_size, block_size, write_mutex);
  }

  GetOutput() = C_;
  return true;
}

void MatmulDoubleSTLTask::Worker(size_t start_step, size_t end_step, size_t grid_size, size_t block_size,
                                 std::mutex &write_mutex) {
  // Обрабатываем диапазон итераций [start_step, end_step)
  for (size_t step_i_j = start_step; step_i_j < end_step; ++step_i_j) {
    const size_t step = step_i_j / (grid_size * grid_size);
    const size_t i = (step_i_j % (grid_size * grid_size)) / grid_size;
    const size_t j = step_i_j % grid_size;

    const size_t root = (i + step) % grid_size;

    std::vector<double> local_block(block_size * block_size, 0.0);

    for (size_t bi = 0; bi < block_size; ++bi) {
      for (size_t bj = 0; bj < block_size; ++bj) {
        double sum = 0.0;
        for (size_t bk = 0; bk < block_size; ++bk) {
          const size_t idx_a = ((i * block_size + bi) * n_) + (root * block_size + bk);
          const size_t idx_b = ((root * block_size + bk) * n_) + (j * block_size + bj);
          sum += A_[idx_a] * B_[idx_b];
        }
        local_block[(bi * block_size) + bj] += sum;
      }
    }

    {
      std::scoped_lock<std::mutex> lock(write_mutex);
      for (size_t bi = 0; bi < block_size; ++bi) {
        for (size_t bj = 0; bj < block_size; ++bj) {
          const size_t idx_c = ((i * block_size + bi) * n_) + (j * block_size + bj);
          C_[idx_c] += local_block[(bi * block_size) + bj];
        }
      }
    }
  }
}

bool MatmulDoubleSTLTask::RunSimpleMultiply() {
  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  const size_t num_threads = std::thread::hardware_concurrency();

  if (n >= num_threads) {
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    const size_t rows_per_thread = n / num_threads;

    for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      const size_t start_row = thread_idx * rows_per_thread;
      const size_t end_row = (thread_idx == num_threads - 1) ? n : start_row + rows_per_thread;

      threads.emplace_back(SimpleMultiplyThread, std::cref(a), std::cref(b), std::ref(c), n, start_row, end_row);
    }

    for (auto &thread : threads) {
      thread.join();
    }
  } else {
    // Простое однопоточное умножение
    SimpleMultiplyThread(a, b, c, n, 0, n);
  }

  return true;
}

bool MatmulDoubleSTLTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_stl
