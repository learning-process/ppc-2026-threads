#include "baranov_a_mult_matrix_fox_algorithm/all/include/ops_all.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "oneapi/tbb.h"

namespace baranov_a_mult_matrix_fox_algorithm_all {

namespace {

void MultiplyBlock(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                   std::vector<double> &output, size_t n, size_t i_start, size_t i_end, size_t j_start, size_t j_end,
                   size_t k_start, size_t k_end) {
  for (size_t i = i_start; i < i_end; ++i) {
    for (size_t j = j_start; j < j_end; ++j) {
      double sum = 0.0;
      for (size_t k = k_start; k < k_end; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] += sum;
    }
  }
}

void MultiplySEQ(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] = sum;
    }
  }
}

#ifdef _OPENMP
void MultiplyOMP(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n) {
#  pragma omp parallel for  // NOLINT(openmp-use-default-none)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] = sum;
    }
  }
}
#endif

#ifdef TBB
void MultiplyTBB(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n) {
  tbb::parallel_for(static_cast<size_t>(0), n, [&](size_t i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] = sum;
    }
  });
}
#endif

void FoxBlockSEQ(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n, size_t block_size) {
  size_t num_blocks = (n + block_size - 1) / block_size;

  // NOLINTNEXTLINE(modernize-use-ranges)
  std::fill(output.begin(), output.end(), 0.0);

  for (size_t bk = 0; bk < num_blocks; ++bk) {
    for (size_t bi = 0; bi < num_blocks; ++bi) {
      for (size_t bj = 0; bj < num_blocks; ++bj) {
        size_t broadcast_block = (bi + bk) % num_blocks;
        size_t i_start = bi * block_size;
        size_t i_end = std::min(i_start + block_size, n);
        size_t j_start = bj * block_size;
        size_t j_end = std::min(j_start + block_size, n);
        size_t k_start = broadcast_block * block_size;
        size_t k_end = std::min(k_start + block_size, n);

        MultiplyBlock(matrix_a, matrix_b, output, n, i_start, i_end, j_start, j_end, k_start, k_end);
      }
    }
  }
}

#ifdef _OPENMP
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void FoxBlockOMP(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n, size_t block_size) {
  size_t num_blocks = (n + block_size - 1) / block_size;

#  pragma omp parallel for  // NOLINT(openmp-use-default-none)
  for (size_t idx = 0; idx < n * n; ++idx) {
    output[idx] = 0.0;
  }

  for (size_t bk = 0; bk < num_blocks; ++bk) {
#  pragma omp parallel for  // NOLINT(openmp-use-default-none)
    for (size_t linear_idx = 0; linear_idx < num_blocks * num_blocks; ++linear_idx) {
      size_t bi = linear_idx / num_blocks;
      size_t bj = linear_idx % num_blocks;
      size_t broadcast_block = (bi + bk) % num_blocks;
      size_t i_start = bi * block_size;
      size_t i_end = std::min(i_start + block_size, n);
      size_t j_start = bj * block_size;
      size_t j_end = std::min(j_start + block_size, n);
      size_t k_start = broadcast_block * block_size;
      size_t k_end = std::min(k_start + block_size, n);

      for (size_t i = i_start; i < i_end; ++i) {
        for (size_t j = j_start; j < j_end; ++j) {
          double sum = 0.0;
          for (size_t k = k_start; k < k_end; ++k) {
            sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
          }
          output[(i * n) + j] += sum;
        }
      }
    }
  }
}
#endif

#ifdef TBB
void FoxBlockTBB(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                 size_t n, size_t block_size) {
  size_t num_blocks = (n + block_size - 1) / block_size;

  tbb::parallel_for(static_cast<size_t>(0), n * n, [&](size_t idx) { output[idx] = 0.0; });

  for (size_t bk = 0; bk < num_blocks; ++bk) {
    tbb::parallel_for(static_cast<size_t>(0), num_blocks * num_blocks, [&](size_t linear_idx) {
      size_t bi = linear_idx / num_blocks;
      size_t bj = linear_idx % num_blocks;
      size_t broadcast_block = (bi + bk) % num_blocks;
      size_t i_start = bi * block_size;
      size_t i_end = std::min(i_start + block_size, n);
      size_t j_start = bj * block_size;
      size_t j_end = std::min(j_start + block_size, n);
      size_t k_start = broadcast_block * block_size;
      size_t k_end = std::min(k_start + block_size, n);

      for (size_t i = i_start; i < i_end; ++i) {
        for (size_t j = j_start; j < j_end; ++j) {
          double sum = 0.0;
          for (size_t k = k_start; k < k_end; ++k) {
            sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
          }
          output[(i * n) + j] += sum;
        }
      }
    });
  }
}
#endif

void MultiplyDispatch(bool use_parallel, const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                      std::vector<double> &output, size_t n) {
  if (!use_parallel) {
    MultiplySEQ(matrix_a, matrix_b, output, n);
    return;
  }

#ifdef TBB
  MultiplyTBB(matrix_a, matrix_b, output, n);
#elifdef _OPENMP  // NOLINT(readability-use-concise-preprocessor-directives)
  MultiplyOMP(matrix_a, matrix_b, output, n);
#else
  MultiplySEQ(matrix_a, matrix_b, output, n);
#endif
}

void FoxBlockDispatch(bool use_parallel, const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                      std::vector<double> &output, size_t n, size_t block_size) {
  if (!use_parallel) {
    FoxBlockSEQ(matrix_a, matrix_b, output, n, block_size);
    return;
  }

#ifdef TBB
  FoxBlockTBB(matrix_a, matrix_b, output, n, block_size);
#elifdef _OPENMP  // NOLINT(readability-use-concise-preprocessor-directives)
  FoxBlockOMP(matrix_a, matrix_b, output, n, block_size);
#else
  FoxBlockSEQ(matrix_a, matrix_b, output, n, block_size);
#endif
}

}  // namespace

BaranovAMultMatrixFoxAlgorithmALL::BaranovAMultMatrixFoxAlgorithmALL(
    const baranov_a_mult_matrix_fox_algorithm::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool BaranovAMultMatrixFoxAlgorithmALL::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool BaranovAMultMatrixFoxAlgorithmALL::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void BaranovAMultMatrixFoxAlgorithmALL::StandardMultiplication(size_t n) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();
  MultiplyDispatch(true, matrix_a, matrix_b, output, n);
}

void BaranovAMultMatrixFoxAlgorithmALL::FoxBlockMultiplication(size_t n, size_t block_size) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();
  FoxBlockDispatch(true, matrix_a, matrix_b, output, n, block_size);
}

bool BaranovAMultMatrixFoxAlgorithmALL::RunImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  size_t n = matrix_size;
  size_t block_size = 64;
  if (n < block_size) {
    StandardMultiplication(n);
  } else {
    FoxBlockMultiplication(n, block_size);
  }
  return true;
}

bool BaranovAMultMatrixFoxAlgorithmALL::PostProcessingImpl() {
  return true;
}

}  // namespace baranov_a_mult_matrix_fox_algorithm_all
