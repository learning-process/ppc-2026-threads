#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/all/include/ops_all.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::RemizovKDenseMatrixMultiplicationCannonAlgorithmAll(
    const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::ValidationImpl() {
  const auto &input_data = GetInput();

  int block_dim = std::get<0>(input_data);
  const auto &mat_a = std::get<1>(input_data);
  const auto &mat_b = std::get<2>(input_data);

  if (block_dim <= 0) {
    return false;
  }
  if (mat_a.empty() || mat_b.empty()) {
    return false;
  }

  size_t n = mat_a.size();
  if (n != mat_a[0].size()) {
    return false;
  }
  if (n != mat_b.size() || n != mat_b[0].size()) {
    return false;
  }

  return (n % static_cast<size_t>(block_dim) == 0);
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::MultiplyBlock(const std::vector<std::vector<double>> &a,
                                                                        const std::vector<std::vector<double>> &b,
                                                                        std::vector<std::vector<double>> &c,
                                                                        int block_size) {
#ifdef _OPENMP
#  pragma omp parallel for collapse(2) schedule(static) default(none) shared(a, b, c, block_size)
#endif
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double acc = 0.0;
      for (int k = 0; k < block_size; ++k) {
        acc += a[i][k] * b[k][j];
      }
      c[i][j] += acc;
    }
  }
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::ShiftBlocksLeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks, int block_count) {
  const unsigned int num_threads = std::max(1U, std::thread::hardware_concurrency());
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  const int rows_per_thread = (block_count + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const int start = static_cast<int>(thread_idx) * rows_per_thread;
    const int end = std::min(start + rows_per_thread, block_count);
    if (start >= end) {
      break;
    }
    threads.emplace_back([&matrix_blocks, block_count, start, end]() {
      for (int i = start; i < end; ++i) {
        auto first = std::move(matrix_blocks[i][0]);
        for (int j = 1; j < block_count; ++j) {
          matrix_blocks[i][j - 1] = std::move(matrix_blocks[i][j]);
        }
        matrix_blocks[i][block_count - 1] = std::move(first);
      }
    });
  }
  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::ShiftBlocksUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks, int block_count) {
  const unsigned int num_threads = std::max(1U, std::thread::hardware_concurrency());
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  const int cols_per_thread = (block_count + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const int start = static_cast<int>(thread_idx) * cols_per_thread;
    const int end = std::min(start + cols_per_thread, block_count);
    if (start >= end) {
      break;
    }
    threads.emplace_back([&matrix_blocks, block_count, start, end]() {
      for (int j = start; j < end; ++j) {
        auto first = std::move(matrix_blocks[0][j]);
        for (int i = 1; i < block_count; ++i) {
          matrix_blocks[i - 1][j] = std::move(matrix_blocks[i][j]);
        }
        matrix_blocks[block_count - 1][j] = std::move(first);
      }
    });
  }
  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::RunCannonCycle(
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks, int block_size, int block_count) {
  for (int step = 0; step < block_count; ++step) {
    tbb::parallel_for(tbb::blocked_range2d<int>(0, block_count, 0, block_count),
                      [&](const tbb::blocked_range2d<int> &r) {
      for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
        for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
          MultiplyBlock(a_blocks[i][j], b_blocks[i][j], c_blocks[i][j], block_size);
        }
      }
    });

    if (step < block_count - 1) {
      ShiftBlocksLeft(a_blocks, block_count);
      ShiftBlocksUp(b_blocks, block_count);
    }
  }
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::InitializeBlocks(
    const std::vector<std::vector<double>> &matrix_a, const std::vector<std::vector<double>> &matrix_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks, int block_size, int block_count) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_count, 0, block_count), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
        const int shift = (i + j) % block_count;
        for (int bi = 0; bi < block_size; ++bi) {
          for (int bj = 0; bj < block_size; ++bj) {
            a_blocks[i][j][bi][bj] = matrix_a[(i * block_size) + bi][(shift * block_size) + bj];
            b_blocks[i][j][bi][bj] = matrix_b[(shift * block_size) + bi][(j * block_size) + bj];
          }
        }
      }
    }
  });
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::AssembleOutput(
    std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks, std::vector<std::vector<double>> &output,
    int block_size, int block_count) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_count, 0, block_count), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
        for (int bi = 0; bi < block_size; ++bi) {
          for (int bj = 0; bj < block_size; ++bj) {
            output[(i * block_size) + bi][(j * block_size) + bj] = c_blocks[i][j][bi][bj];
          }
        }
      }
    }
  });
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::RunImpl() {
  const auto &params = GetInput();

  const int block_dim = std::get<0>(params);
  const auto &source_a = std::get<1>(params);
  const auto &source_b = std::get<2>(params);

  const int matrix_size = static_cast<int>(source_a.size());
  const int blocks_per_dim = matrix_size / block_dim;

  using Block4D = std::vector<std::vector<std::vector<std::vector<double>>>>;
  Block4D blocks_a(blocks_per_dim, std::vector<std::vector<std::vector<double>>>(
                                       blocks_per_dim, std::vector<std::vector<double>>(
                                                           block_dim, std::vector<double>(block_dim, 0.0))));
  Block4D blocks_b = blocks_a;
  Block4D blocks_c = blocks_a;

  InitializeBlocks(source_a, source_b, blocks_a, blocks_b, block_dim, blocks_per_dim);
  RunCannonCycle(blocks_a, blocks_b, blocks_c, block_dim, blocks_per_dim);

  std::vector<std::vector<double>> result(matrix_size, std::vector<double>(matrix_size, 0.0));
  AssembleOutput(blocks_c, result, block_dim, blocks_per_dim);

  GetOutput() = std::move(result);
  return true;
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmAll::PostProcessingImpl() {
  return true;
}

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
