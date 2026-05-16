#include "sinev_a_mult_matrix_fox_algorithm/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmTBB::SinevAMultMatrixFoxAlgorithmTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmTBB::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmTBB::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

// Добавляем static к определениям
void SinevAMultMatrixFoxAlgorithmTBB::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, n, 0, n), [&](const tbb::blocked_range2d<size_t> &r) {
    for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
      for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int bi = r.rows().begin(); bi < r.rows().end(); ++bi) {
      for (int bj = r.cols().begin(); bj < r.cols().end(); ++bj) {
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
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int bi = r.rows().begin(); bi < r.rows().end(); ++bi) {
      for (int bj = r.cols().begin(); bj < r.cols().end(); ++bj) {
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
  });
}

void SinevAMultMatrixFoxAlgorithmTBB::MultiplyBlocks(const std::vector<double> &blocks_a,
                                                     const std::vector<double> &blocks_b, std::vector<double> &blocks_c,
                                                     size_t bs, size_t a_off, size_t b_off, size_t c_off) {
  for (size_t ii = 0; ii < bs; ++ii) {
    for (size_t kk = 0; kk < bs; ++kk) {
      const double val = blocks_a[a_off + (ii * bs) + kk];
      const size_t b_base = b_off + (kk * bs);
      const size_t c_base = c_off + (ii * bs);
      for (size_t jj = 0; jj < bs; ++jj) {
        blocks_c[c_base + jj] += val * blocks_b[b_base + jj];
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmTBB::FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                                              std::vector<double> &blocks_c, size_t bs, int q, int step) {
  const size_t block_size = bs * bs;

  tbb::parallel_for(tbb::blocked_range2d<int>(0, q, 0, q), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
        const int k = (i + step) % q;

        const size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size;
        const size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size;
        const size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size;

        MultiplyBlocks(blocks_a, blocks_b, blocks_c, bs, a_off, b_off, c_off);
      }
    }
  });
}

bool SinevAMultMatrixFoxAlgorithmTBB::RunImpl() {
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

  size_t bs = ChooseBlockSize(n);
  const int actual_q = static_cast<int>(n / bs);

  const auto total_blocks = static_cast<size_t>(actual_q) * static_cast<size_t>(actual_q);
  const auto block_elements = bs * bs;

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

size_t SinevAMultMatrixFoxAlgorithmTBB::ChooseBlockSize(size_t n) {
  size_t bs = 1;
  const auto sqrt_n = static_cast<size_t>(std::sqrt(static_cast<double>(n)));
  for (size_t div = sqrt_n; div >= 1; --div) {
    if (n % div == 0) {
      bs = div;
      break;
    }
  }
  return bs;
}

bool SinevAMultMatrixFoxAlgorithmTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
