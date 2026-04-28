#include "baranov_a_mult_matrix_fox_algorithm/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "oneapi/tbb.h"

namespace baranov_a_mult_matrix_fox_algorithm_tbb {

BaranovAMultMatrixFoxAlgorithmTBB::BaranovAMultMatrixFoxAlgorithmTBB(
    const baranov_a_mult_matrix_fox_algorithm::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool BaranovAMultMatrixFoxAlgorithmTBB::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool BaranovAMultMatrixFoxAlgorithmTBB::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void BaranovAMultMatrixFoxAlgorithmTBB::StandardMultiplication(size_t n) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();

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

void BaranovAMultMatrixFoxAlgorithmTBB::FoxBlockMultiplication(size_t n, size_t block_size) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();

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

bool BaranovAMultMatrixFoxAlgorithmTBB::RunImpl() {
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

bool BaranovAMultMatrixFoxAlgorithmTBB::PostProcessingImpl() {
  return true;
}

}  // namespace baranov_a_mult_matrix_fox_algorithm_tbb
