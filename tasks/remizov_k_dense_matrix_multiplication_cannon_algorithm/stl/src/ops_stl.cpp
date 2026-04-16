#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>
#include <functional>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

template <typename IndexType, typename Func>
static void ParallelFor(IndexType begin, IndexType end, Func&& func) {
  const std::size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
  const IndexType range_length = end - begin;
  if (range_length <= 0) return;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  IndexType chunk_size = (range_length + num_threads - 1) / num_threads;
  IndexType start = begin;

  for (std::size_t t = 0; t < num_threads; ++t) {
    IndexType chunk_end = std::min(end, start + chunk_size);
    if (start >= chunk_end) break;

    threads.emplace_back([start, chunk_end, &func]() {
      for (IndexType i = start; i < chunk_end; ++i) {
        func(i);
      }
    });
    start = chunk_end;
  }

  for (auto& th : threads) {
    if (th.joinable()) th.join();
  }
}

template <typename Func>
static void ParallelFor2D(int rows_begin, int rows_end, int cols_begin, int cols_end, Func&& func) {
  const int rows = rows_end - rows_begin;
  const int cols = cols_end - cols_begin;
  const int total = rows * cols;
  if (total <= 0) return;

  ParallelFor(0, total, [&](int linear_idx) {
    int i = rows_begin + linear_idx / cols;
    int j = cols_begin + linear_idx % cols;
    func(i, j);
  });
}


RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::RemizovKDenseMatrixMultiplicationCannonAlgorithmStl(
    const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}


bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::ValidationImpl() {
  const auto &input_data = GetInput();
  int block_dim = std::get<0>(input_data);
  const auto &mat_a = std::get<1>(input_data);
  const auto &mat_b = std::get<2>(input_data);

  if (block_dim <= 0) return false;
  if (mat_a.empty() || mat_b.empty()) return false;

  size_t n = mat_a.size();
  if (n != mat_a[0].size()) return false;
  if (n != mat_b.size() || n != mat_b[0].size()) return false;

  return (n % static_cast<size_t>(block_dim) == 0);
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}


void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::MultiplyBlock(
    const std::vector<std::vector<double>> &a,
    const std::vector<std::vector<double>> &b,
    std::vector<std::vector<double>> &c,
    int block_size) {
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

void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::ShiftBlocksLeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks,
    int block_count) {
  ParallelFor(0, block_count, [&](int i) {
    auto first = std::move(matrix_blocks[i][0]);
    for (int j = 1; j < block_count; ++j) {
      matrix_blocks[i][j - 1] = std::move(matrix_blocks[i][j]);
    }
    matrix_blocks[i][block_count - 1] = std::move(first);
  });
}


void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::ShiftBlocksUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks,
    int block_count) {
  ParallelFor(0, block_count, [&](int j) {
    auto first = std::move(matrix_blocks[0][j]);
    for (int i = 1; i < block_count; ++i) {
      matrix_blocks[i - 1][j] = std::move(matrix_blocks[i][j]);
    }
    matrix_blocks[block_count - 1][j] = std::move(first);
  });
}


void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::RunCannonCycle(
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks,
    int block_size,
    int block_count) {
  for (int step = 0; step < block_count; ++step) {
    ParallelFor2D(0, block_count, 0, block_count,
                  [&](int i, int j) {
                    MultiplyBlock(a_blocks[i][j], b_blocks[i][j], c_blocks[i][j], block_size);
                  });

    if (step < block_count - 1) {
      ShiftBlocksLeft(a_blocks, block_count);
      ShiftBlocksUp(b_blocks, block_count);
    }
  }
}


void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::InitializeBlocks(
    const std::vector<std::vector<double>> &matrix_a,
    const std::vector<std::vector<double>> &matrix_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks,
    int block_size,
    int block_count) {
  ParallelFor2D(0, block_count, 0, block_count,
                [&](int i, int j) {
                  int shift = (i + j) % block_count;
                  for (int bi = 0; bi < block_size; ++bi) {
                    for (int bj = 0; bj < block_size; ++bj) {
                      a_blocks[i][j][bi][bj] = matrix_a[(i * block_size) + bi][(shift * block_size) + bj];
                      b_blocks[i][j][bi][bj] = matrix_b[(shift * block_size) + bi][(j * block_size) + bj];
                    }
                  }
                });
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::AssembleOutput(
    std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks,
    std::vector<std::vector<double>> &output,
    int block_size,
    int block_count) {
  ParallelFor2D(0, block_count, 0, block_count,
                [&](int i, int j) {
                  for (int bi = 0; bi < block_size; ++bi) {
                    for (int bj = 0; bj < block_size; ++bj) {
                      output[(i * block_size) + bi][(j * block_size) + bj] = c_blocks[i][j][bi][bj];
                    }
                  }
                });
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::RunImpl() {
  const auto &params = GetInput();
  int block_dim = std::get<0>(params);
  const auto &source_a = std::get<1>(params);
  const auto &source_b = std::get<2>(params);

  int matrix_size = static_cast<int>(source_a.size());
  int blocks_per_dim = matrix_size / block_dim;

  using Block4D = std::vector<std::vector<std::vector<std::vector<double>>>>;
  Block4D blocks_a(blocks_per_dim,
                   std::vector<std::vector<std::vector<double>>>(
                       blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));
  Block4D blocks_b(blocks_per_dim,
                   std::vector<std::vector<std::vector<double>>>(
                       blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));
  Block4D blocks_c(blocks_per_dim,
                   std::vector<std::vector<std::vector<double>>>(
                       blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));

  InitializeBlocks(source_a, source_b, blocks_a, blocks_b, block_dim, blocks_per_dim);
  RunCannonCycle(blocks_a, blocks_b, blocks_c, block_dim, blocks_per_dim);

  std::vector<std::vector<double>> result(matrix_size, std::vector<double>(matrix_size, 0.0));
  AssembleOutput(blocks_c, result, block_dim, blocks_per_dim);

  GetOutput() = std::move(result);
  return true;
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::PostProcessingImpl() {
  return true;
}

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
