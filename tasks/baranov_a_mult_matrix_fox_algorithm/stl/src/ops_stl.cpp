#include "baranov_a_mult_matrix_fox_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_stl {

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

void MultiplyRowRange(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                      std::vector<double> &output, size_t n, size_t start_i, size_t end_i) {
  for (size_t i = start_i; i < end_i; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < n; ++k) {
        sum += matrix_a[(i * n) + k] * matrix_b[(k * n) + j];
      }
      output[(i * n) + j] = sum;
    }
  }
}

void ParallelRowMultiplication(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                               std::vector<double> &output, size_t n, size_t start_i, size_t end_i) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }

  std::vector<std::thread> threads;
  size_t chunk_size = (end_i - start_i + num_threads - 1) / num_threads;

  for (unsigned int tid = 0; tid < num_threads; ++tid) {
    size_t i_start = start_i + (tid * chunk_size);
    size_t i_end_local = std::min(i_start + chunk_size, end_i);
    if (i_start >= end_i) {
      break;
    }

    threads.emplace_back(
        [&, i_start, i_end_local]() { MultiplyRowRange(matrix_a, matrix_b, output, n, i_start, i_end_local); });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void ProcessBlockRange(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                       std::vector<double> &output, size_t n, size_t block_k, size_t num_blocks, size_t block_size,
                       const std::vector<size_t> &block_indices, size_t start_idx, size_t end_idx) {
  for (size_t idx = start_idx; idx < end_idx; ++idx) {
    size_t linear_idx = block_indices[idx];
    size_t block_i = linear_idx / num_blocks;
    size_t block_j = linear_idx % num_blocks;

    size_t broadcast_block = (block_i + block_k) % num_blocks;

    size_t i_start = block_i * block_size;
    size_t i_end = std::min(i_start + block_size, n);
    size_t j_start = block_j * block_size;
    size_t j_end = std::min(j_start + block_size, n);
    size_t k_start = broadcast_block * block_size;
    size_t k_end = std::min(k_start + block_size, n);

    MultiplyBlock(matrix_a, matrix_b, output, n, i_start, i_end, j_start, j_end, k_start, k_end);
  }
}

void ParallelBlockProcessing(const std::vector<double> &matrix_a, const std::vector<double> &matrix_b,
                             std::vector<double> &output, size_t n, size_t block_k, size_t num_blocks,
                             size_t block_size, const std::vector<size_t> &block_indices) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }

  std::vector<std::thread> threads;
  size_t chunk_size = (block_indices.size() + num_threads - 1) / num_threads;

  for (unsigned int thread_id = 0; thread_id < num_threads; ++thread_id) {
    size_t start_idx = thread_id * chunk_size;
    size_t end_idx = std::min(start_idx + chunk_size, block_indices.size());
    if (start_idx >= block_indices.size()) {
      break;
    }

    threads.emplace_back([&, start_idx, end_idx]() {
      ProcessBlockRange(matrix_a, matrix_b, output, n, block_k, num_blocks, block_size, block_indices, start_idx,
                        end_idx);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

}  // namespace

BaranovAMultMatrixFoxAlgorithmSTL::BaranovAMultMatrixFoxAlgorithmSTL(
    const baranov_a_mult_matrix_fox_algorithm::InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool BaranovAMultMatrixFoxAlgorithmSTL::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool BaranovAMultMatrixFoxAlgorithmSTL::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void BaranovAMultMatrixFoxAlgorithmSTL::StandardMultiplication(size_t n) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();
  ParallelRowMultiplication(matrix_a, matrix_b, output, n, 0, n);
}

void BaranovAMultMatrixFoxAlgorithmSTL::FoxBlockMultiplication(size_t n, size_t block_size) {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  auto &output = GetOutput();

  size_t num_blocks = (n + block_size - 1) / block_size;

  std::ranges::fill(output, 0.0);

  std::vector<size_t> block_indices(num_blocks * num_blocks);
  for (size_t idx = 0; idx < num_blocks * num_blocks; ++idx) {
    block_indices[idx] = idx;
  }

  for (size_t block_k = 0; block_k < num_blocks; ++block_k) {
    ParallelBlockProcessing(matrix_a, matrix_b, output, n, block_k, num_blocks, block_size, block_indices);
  }
}

bool BaranovAMultMatrixFoxAlgorithmSTL::RunImpl() {
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

bool BaranovAMultMatrixFoxAlgorithmSTL::PostProcessingImpl() {
  return true;
}

}  // namespace baranov_a_mult_matrix_fox_algorithm_stl
