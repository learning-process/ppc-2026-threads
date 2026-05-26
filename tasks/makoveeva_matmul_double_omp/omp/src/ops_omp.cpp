#include "makoveeva_matmul_double_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "makoveeva_matmul_double_omp/common/include/common.hpp"

namespace makoveeva_matmul_double_omp {

namespace {

// Выбирает размер блока в зависимости от размера матрицы
[[nodiscard]] size_t SelectBlockSize(size_t n) {
  // Используем степени двойки для хорошей локальности кэша
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

// Декодирует одномерный индекс в трёхмерный индекс (step, i, j)
void DecodeIndex(size_t step_i_j, size_t grid_size, size_t &step, size_t &i, size_t &j) {
  step = step_i_j / (grid_size * grid_size);
  i = (step_i_j % (grid_size * grid_size)) / grid_size;
  j = step_i_j % grid_size;
}

// Вычисляет root блок для алгоритма Фокса
[[nodiscard]] size_t ComputeRoot(size_t i, size_t step, size_t grid_size) {
  return (i + step) % grid_size;
}

// Умножает блок A[i][root] на блок B[root][j] и сохраняет в local_block
void MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &local_block,
                    size_t i, size_t root, size_t j, size_t block_size, size_t n) {
  for (size_t bi = 0; bi < block_size; ++bi) {
    for (size_t bj = 0; bj < block_size; ++bj) {
      double sum = 0.0;
      for (size_t bk = 0; bk < block_size; ++bk) {
        const size_t idx_a = ((i * block_size + bi) * n) + (root * block_size + bk);
        const size_t idx_b = ((root * block_size + bk) * n) + (j * block_size + bj);
        sum += a[idx_a] * b[idx_b];
      }
      local_block[(bi * block_size) + bj] += sum;
    }
  }
}

// Добавляет результат из local_block в матрицу C
void AddBlockToResult(std::vector<double> &c, const std::vector<double> &local_block, size_t i, size_t j,
                      size_t block_size, size_t n) {
  for (size_t bi = 0; bi < block_size; ++bi) {
    for (size_t bj = 0; bj < block_size; ++bj) {
      const size_t idx_c = ((i * block_size + bi) * n) + (j * block_size + bj);
      c[idx_c] += local_block[(bi * block_size) + bj];
    }
  }
}

}  // namespace

MatmulDoubleOMPTask::MatmulDoubleOMPTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleOMPTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleOMPTask::PreProcessingImpl() {
  const auto &input = GetInput();
  n_ = std::get<0>(input);
  a_ = std::get<1>(input);
  b_ = std::get<2>(input);
  c_.assign(n_ * n_, 0.0);

  return true;
}

bool MatmulDoubleOMPTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  const size_t n = n_;
  const auto &a = a_;
  const auto &b = b_;
  auto &c = c_;

  // Выбираем размер блока для оптимальной локальности кэша
  const size_t block_size = SelectBlockSize(n);

  // Проверяем что матрица делится нацело на размер блока
  if (n % block_size != 0) {
    return RunSimpleMultiply();
  }

  const size_t grid_size = n / block_size;

  // Алгоритм Фокса: все итерации параллелизируются в один цикл
#pragma omp parallel for default(none) shared(a, b, c, n, block_size, grid_size)
  for (size_t step_i_j = 0; step_i_j < grid_size * grid_size * grid_size; ++step_i_j) {
    size_t step = 0;
    size_t i = 0;
    size_t j = 0;
    DecodeIndex(step_i_j, grid_size, step, i, j);

    // Источник блока A: диагональный сдвиг на step позиций
    const size_t root = ComputeRoot(i, step, grid_size);

    // Локальный буфер для накопления результатов блока C[i][j]
    std::vector<double> local_block(block_size * block_size, 0.0);

    // Умножение блока A[i][root] на блок B[root][j]
    MultiplyBlocks(a, b, local_block, i, root, j, block_size, n);

    // Безопасно добавляем результат в матрицу C
#pragma omp critical
    {
      AddBlockToResult(c, local_block, i, j, block_size, n);
    }
  }

  GetOutput() = c_;
  return true;
}

bool MatmulDoubleOMPTask::RunSimpleMultiply() {
  const size_t n = n_;
  const auto &a = a_;
  const auto &b = b_;
  auto &c = c_;

#pragma omp parallel for collapse(2) default(none) shared(a, b, c, n)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] = sum;
    }
  }

  return true;
}

bool MatmulDoubleOMPTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_omp
