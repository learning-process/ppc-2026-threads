#include "sinev_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include <omp.h>

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmOMP::SinevAMultMatrixFoxAlgorithmOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmOMP::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmOMP::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

// Простое умножение для маленьких матриц
void SinevAMultMatrixFoxAlgorithmOMP::SimpleMultiply(size_t n, const std::vector<double>& A, 
                    const std::vector<double>& B, std::vector<double>& C) {
  #pragma omp parallel for collapse(2)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += A[i * n + k] * B[k * n + j];
      }
      C[i * n + j] = sum;
    }
  }
}

bool SinevAMultMatrixFoxAlgorithmOMP::RunImpl() {
  const auto &input = GetInput();
  size_t n = std::get<0>(input);
  const auto &A = std::get<1>(input);
  const auto &B = std::get<2>(input);
  auto &C = GetOutput();

  // Для маленьких матриц используем простое умножение
  if (n <= 8) {
    SimpleMultiply(n, A, B, C);
    return true;
  }

  int num_threads = omp_get_max_threads();
  int q = static_cast<int>(std::sqrt(num_threads));
  while (q * q > num_threads) q--;
  if (q < 1) q = 1;

  size_t bs = 1;
  for (size_t div = static_cast<size_t>(std::sqrt(n)); div >= 1; --div) {
    if (n % div == 0) {
      bs = div;
      break;
    }
  }
  
  // Пересчитываем q под выбранный размер блока
  q = n / bs;

  std::vector<double> blocksA(q * q * bs * bs);
  std::vector<double> blocksB(q * q * bs * bs);
  std::vector<double> blocksC(q * q * bs * bs, 0.0);

  // разложение матрицы на блоки
  #pragma omp parallel for collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          blocksA[block_off + i * bs + j] = A[((bi * bs + i) * n) + (bj * bs + j)];
          blocksB[block_off + i * bs + j] = B[((bi * bs + i) * n) + (bj * bs + j)];
        }
      }
    }
  }

  // алгоритм фокса
  for (int step = 0; step < q; ++step) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < q; ++i) {
      for (int j = 0; j < q; ++j) {
        int k = (i + step) % q; 
        
        int a_off = ((i * q) + k) * (bs * bs);
        int b_off = ((k * q) + j) * (bs * bs);
        int c_off = ((i * q) + j) * (bs * bs);
        
        for (size_t ii = 0; ii < bs; ++ii) {
          for (size_t kk = 0; kk < bs; ++kk) {
            double val = blocksA[a_off + ii * bs + kk];
            for (size_t jj = 0; jj < bs; ++jj) {
              blocksC[c_off + ii * bs + jj] += val * blocksB[b_off + kk * bs + jj];
            }
          }
        }
      }
    }
  }

  // сбор результата
  #pragma omp parallel for collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      int block_off = ((bi * q) + bj) * (bs * bs);
      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          C[((bi * bs + i) * n) + (bj * bs + j)] = blocksC[block_off + i * bs + j];
        }
      }
    }
  }

  return true;
}

bool SinevAMultMatrixFoxAlgorithmOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
