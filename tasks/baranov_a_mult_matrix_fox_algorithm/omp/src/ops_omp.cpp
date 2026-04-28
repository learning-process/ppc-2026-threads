#include "baranov_a_mult_matrix_fox_algorithm/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace {

void ProcessBlock(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b, std::vector<double> &output,
                  size_t n, size_t i_start, size_t i_end, size_t j_start, size_t j_end, size_t k_start, size_t k_end) {
  for (size_t i = i_start; i < i_end; ++i) {
    for (size_t j = j_start; j < j_end; ++j) {
      double sum = 0.0;
      for (size_t k = k_start; k < k_end; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
#pragma omp atomic
      output[(i * n) + j] += sum;
    }
  }
}

}  // namespace

namespace baranov_a_mult_matrix_fox_algorithm_omp {

BaranovAMultMatrixFoxAlgorithmOMP::BaranovAMultMatrixFoxAlgorithmOMP(
    const baranov_a_mult_matrix_fox_algorithm::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool BaranovAMultMatrixFoxAlgorithmOMP::ValidationImpl() {
  const auto &input = GetInput();
  size_t matrix_size = std::get<0>(input);
  const auto &matrix_a = std::get<1>(input);
  const auto &matrix_b = std::get<2>(input);

  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool BaranovAMultMatrixFoxAlgorithmOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  size_t matrix_size = std::get<0>(input);

  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void BaranovAMultMatrixFoxAlgorithmOMP::StandardMultiplication(size_t n) {
  const auto &input = GetInput();
  const auto &matrix_a = std::get<1>(input);
  const auto &matrix_b = std::get<2>(input);
  auto &output = GetOutput();

#pragma omp parallel for collapse(2) default(none) shared(matrix_a, matrix_b, output, n)
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

void BaranovAMultMatrixFoxAlgorithmOMP::FoxBlockMultiplication(size_t n, size_t block_size) {
  const auto &input = GetInput();
  const auto &matrix_a = std::get<1>(input);
  const auto &matrix_b = std::get<2>(input);
  auto &output = GetOutput();

  size_t num_blocks = (n + block_size - 1) / block_size;

#pragma omp parallel for default(none) shared(output, n)
  for (size_t idx = 0; idx < n * n; ++idx) {
    output[idx] = 0.0;
  }

  for (size_t bk = 0; bk < num_blocks; ++bk) {
#pragma omp parallel for collapse(2) default(none) shared(matrix_a, matrix_b, output, n, block_size, num_blocks, bk)
    for (size_t bi = 0; bi < num_blocks; ++bi) {
      for (size_t bj = 0; bj < num_blocks; ++bj) {
        size_t broadcast_block = (bi + bk) % num_blocks;

        size_t i_start = bi * block_size;
        size_t i_end = std::min(i_start + block_size, n);
        size_t j_start = bj * block_size;
        size_t j_end = std::min(j_start + block_size, n);
        size_t k_start = broadcast_block * block_size;
        size_t k_end = std::min(k_start + block_size, n);

        ProcessBlock(matrix_a, matrix_b, output, n, i_start, i_end, j_start, j_end, k_start, k_end);
      }
    }
  }
}

bool BaranovAMultMatrixFoxAlgorithmOMP::RunImpl() {
  const auto &input = GetInput();
  size_t n = std::get<0>(input);

  int num_threads = omp_get_max_threads();
  size_t block_size = 64;

  if (num_threads > 1) {
    block_size = std::max(static_cast<size_t>(32), n / static_cast<size_t>(num_threads * 2));
    block_size = std::min(block_size, static_cast<size_t>(128));
  }

  if (n < block_size) {
    StandardMultiplication(n);
  } else {
    FoxBlockMultiplication(n, block_size);
  }

  return true;
}

bool BaranovAMultMatrixFoxAlgorithmOMP::PostProcessingImpl() {
  return true;
}

}  // namespace baranov_a_mult_matrix_fox_algorithm_omp
