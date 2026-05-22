#include "sinev_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

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

void SinevAMultMatrixFoxAlgorithmOMP::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      c[(i * n) + j] = sum;
    }
  }
}

void SinevAMultMatrixFoxAlgorithmOMP::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_off = (static_cast<size_t>((bi * q) + bj)) * (bs * bs);
      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
          const size_t dst_idx = block_off + (i * bs) + j;
          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmOMP::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_off = (static_cast<size_t>((bi * q) + bj)) * (bs * bs);
      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = block_off + (i * bs) + j;
          const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmOMP::FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                                              std::vector<double> &blocks_c, size_t bs, int q, int step) {
  const size_t block_size_bytes = bs * bs;
#pragma omp parallel for default(none) shared(blocks_a, blocks_b, blocks_c, bs, q, step, block_size_bytes) collapse(2)
  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {
      const int k = (i + step) % q;

      const size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size_bytes;
      const size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size_bytes;
      const size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size_bytes;

      for (size_t ii = 0; ii < bs; ++ii) {
        for (size_t kk = 0; kk < bs; ++kk) {
          const double val = blocks_a[a_off + (ii * bs) + kk];
          for (size_t jj = 0; jj < bs; ++jj) {
            blocks_c[c_off + (ii * bs) + jj] += val * blocks_b[b_off + (kk * bs) + jj];
          }
        }
      }
    }
  }
}

bool SinevAMultMatrixFoxAlgorithmOMP::RunImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);
  auto &c = GetOutput();

  // Для маленьких матриц используем простое умножение
  if (n <= 8) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = 1;
  auto sqrt_n = static_cast<size_t>(std::sqrt(static_cast<double>(n)));
  for (size_t div = sqrt_n; div >= 1; --div) {
    if (n % div == 0) {
      bs = div;
      break;
    }
  }

  const int actual_q = static_cast<int>(n / bs);

  auto total_blocks = static_cast<size_t>(actual_q) * static_cast<size_t>(actual_q);
  auto block_elements = bs * bs;

  std::vector<double> blocks_a(total_blocks * block_elements);
  std::vector<double> blocks_b(total_blocks * block_elements);
  std::vector<double> blocks_c(total_blocks * block_elements, 0.0);

  DecomposeToBlocks(a, blocks_a, n, bs, actual_q);
  DecomposeToBlocks(b, blocks_b, n, bs, actual_q);

  for (int step = 0; step < actual_q; ++step) {
    FoxStep(blocks_a, blocks_b, blocks_c, bs, actual_q, step);
  }

  AssembleFromBlocks(blocks_c, c, n, bs, actual_q);

  return true;
}

bool SinevAMultMatrixFoxAlgorithmOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
