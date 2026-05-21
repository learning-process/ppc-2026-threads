#include "makoveeva_matmul_double_tbb/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/mutex.h>
#include <oneapi/tbb/parallel_for.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_tbb/common/include/common.hpp"

namespace makoveeva_matmul_double_tbb {

namespace {

// Выбирает размер блока в зависимости от размера матрицы
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

// Вычисляет умножение блока A[i][root] на блок B[root][j]
void ComputeBlock(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                  std::vector<double> &local_block, size_t i, size_t j, size_t root,
                  size_t block_size, size_t n) {
  for (size_t bi = 0; bi < block_size; ++bi) {
    for (size_t bj = 0; bj < block_size; ++bj) {
      double sum = 0.0;
      for (size_t bk = 0; bk < block_size; ++bk) {
        const size_t idx_a = ((i * block_size + bi) * n) + (root * block_size + bk);
        const size_t idx_b = ((root * block_size + bk) * n) + (j * block_size + bj);
        sum += matrix_a[idx_a] * matrix_b[idx_b];
      }
      local_block[(bi * block_size) + bj] += sum;
    }
  }
}

// Безопасно добавляет результат из local_block в матрицу C
void AccumulateResult(std::vector<double> &matrix_c, const std::vector<double> &local_block,
                      size_t i, size_t j, size_t block_size, size_t n,
                      oneapi::tbb::mutex &write_mutex) {
  oneapi::tbb::mutex::scoped_lock lock(write_mutex);
  for (size_t bi = 0; bi < block_size; ++bi) {
    for (size_t bj = 0; bj < block_size; ++bj) {
      const size_t idx_c = ((i * block_size + bi) * n) + (j * block_size + bj);
      matrix_c[idx_c] += local_block[(bi * block_size) + bj];
    }
  }
}

}  // namespace

MatmulDoubleTBBTask::MatmulDoubleTBBTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleTBBTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleTBBTask::PreProcessingImpl() {
  const auto &input = GetInput();
  n_ = std::get<0>(input);
  A_ = std::get<1>(input);
  B_ = std::get<2>(input);
  C_.assign(n_ * n_, 0.0);

  return true;
}

bool MatmulDoubleTBBTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  // Выбираем размер блока для оптимальной локальности кэша
  const size_t block_size = SelectBlockSize(n);

  // Проверяем что матрица делится нацело на размер блока
  if (n % block_size != 0) {
    return RunSimpleMultiply();
  }

  const size_t grid_size = n / block_size;

  // Создаём mutex для синхронизации доступа к матрице C
  oneapi::tbb::mutex write_mutex;

  // ========================================================================
  // ОСНОВНОЙ ПАРАЛЛЕЛЬНЫЙ ЦИКЛ - АЛГОРИТМ ФОКСА С TBB
  // ========================================================================
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, grid_size * grid_size * grid_size),
      [&](const oneapi::tbb::blocked_range<size_t> &range) {
        for (size_t step_i_j = range.begin(); step_i_j != range.end(); ++step_i_j) {
          // ЭТАП 1: Декодирование одномерного индекса в трёхмерный (step, i, j)
          const size_t step = step_i_j / (grid_size * grid_size);
          const size_t i = (step_i_j % (grid_size * grid_size)) / grid_size;
          const size_t j = step_i_j % grid_size;

          // ЭТАП 2: Вычисление root блока для алгоритма Фокса
          // root = (i + step) % grid_size - КЛЮЧЕВАЯ ЛИНИЯ!
          const size_t root = (i + step) % grid_size;

          // ЭТАП 3: Создание локального буфера (не общего!)
          std::vector<double> local_block(block_size * block_size, 0.0);

          // ЭТАП 4: Умножение блока A[i][root] на блок B[root][j]
          ComputeBlock(a, b, local_block, i, j, root, block_size, n);

          // ЭТАП 5: Безопасное добавление результата в матрицу C
          AccumulateResult(c, local_block, i, j, block_size, n, write_mutex);
        }
      });

  GetOutput() = C_;
  return true;
}

bool MatmulDoubleTBBTask::RunSimpleMultiply() {
  const size_t n = n_;
  const auto &a = A_;
  const auto &b = B_;
  auto &c = C_;

  // Простое параллельное умножение для маленьких матриц
  oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, n),
      [&](const oneapi::tbb::blocked_range<size_t> &range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
          for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
              sum += a[(i * n) + k] * b[(k * n) + j];
            }
            c[(i * n) + j] = sum;
          }
        }
      });

  return true;
}

bool MatmulDoubleTBBTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_tbb
