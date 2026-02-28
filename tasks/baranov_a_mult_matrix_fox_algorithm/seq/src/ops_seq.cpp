#include "baranov_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace baranov_a_mult_matrix_fox_algorithm_seq {

using namespace baranov_a_mult_matrix_fox_algorithm;

BaranovAMultMatrixFoxAlgorithmSEQ::BaranovAMultMatrixFoxAlgorithmSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool BaranovAMultMatrixFoxAlgorithmSEQ::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool BaranovAMultMatrixFoxAlgorithmSEQ::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void BaranovAMultMatrixFoxAlgorithmSEQ::StandardMultiplication(size_t n) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[i * n + k] * matrix_b[k * n + j];
      }
      output[i * n + j] = sum;
    }
  }
}

void BaranovAMultMatrixFoxAlgorithmSEQ::FoxBlockMultiplication(size_t n, size_t block_size) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();

  size_t num_blocks = n / block_size;
  if (n % block_size != 0) {
    num_blocks++;
  }
  std::fill(output.begin(), output.end(), 0.0);

  for (size_t bi = 0; bi < num_blocks; ++bi) {
    for (size_t bj = 0; bj < num_blocks; ++bj) {
      for (size_t bk = 0; bk < num_blocks; ++bk) {
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
              sum += matrix_a[i * n + k] * matrix_b[k * n + j];
            }
            output[i * n + j] += sum;
          }
        }
      }
    }
  }
}

bool BaranovAMultMatrixFoxAlgorithmSEQ::RunImpl() {
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

bool BaranovAMultMatrixFoxAlgorithmSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace baranov_a_mult_matrix_fox_algorithm_seq
