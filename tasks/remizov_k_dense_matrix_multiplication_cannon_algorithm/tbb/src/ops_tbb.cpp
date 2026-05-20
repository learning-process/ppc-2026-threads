#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb(
    const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::ValidationImpl() {
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

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::MultiplyBlock(const std::vector<std::vector<double>> &a,
                                                                        const std::vector<std::vector<double>> &b,
                                                                        std::vector<std::vector<double>> &c,
                                                                        int block_size) {
  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      double accumulator = 0.0;
      for (int k = 0; k < block_size; ++k) {
        accumulator += a[i][k] * b[k][j];
      }
      c[i][j] += accumulator;
    }
  }
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::ShiftBlocksLeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks, int block_count) {
  // Parallelize over rows; each row shift is independent
  tbb::parallel_for(0, block_count, [&](int i) {
    auto first_element = std::move(matrix_blocks[i][0]);
    for (int j = 1; j < block_count; ++j) {
      matrix_blocks[i][j - 1] = std::move(matrix_blocks[i][j]);
    }
    matrix_blocks[i][block_count - 1] = std::move(first_element);
  });
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::ShiftBlocksUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks, int block_count) {
  // Parallelize over columns; each column shift is independent
  tbb::parallel_for(0, block_count, [&](int j) {
    auto first_element = std::move(matrix_blocks[0][j]);
    for (int i = 1; i < block_count; ++i) {
      matrix_blocks[i - 1][j] = std::move(matrix_blocks[i][j]);
    }
    matrix_blocks[block_count - 1][j] = std::move(first_element);
  });
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::RunCannonCycle(
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks, int block_size, int block_count) {
  for (int step = 0; step < block_count; ++step) {
    // Parallel multiplication of all block pairs
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

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::InitializeBlocks(
    const std::vector<std::vector<double>> &matrix_a, const std::vector<std::vector<double>> &matrix_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
    std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks, int block_size, int block_count) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, block_count, 0, block_count), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
        int shift_value = (i + j) % block_count;
        for (int bi = 0; bi < block_size; ++bi) {
          for (int bj = 0; bj < block_size; ++bj) {
            a_blocks[i][j][bi][bj] = matrix_a[(i * block_size) + bi][(shift_value * block_size) + bj];
            b_blocks[i][j][bi][bj] = matrix_b[(shift_value * block_size) + bi][(j * block_size) + bj];
          }
        }
      }
    }
  });
}

void RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::AssembleOutput(
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

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::RunImpl() {
  const auto &params = GetInput();

  int block_dim = std::get<0>(params);
  const auto &source_a = std::get<1>(params);
  const auto &source_b = std::get<2>(params);

  int matrix_size = static_cast<int>(source_a.size());
  int blocks_per_dim = matrix_size / block_dim;

  std::vector<std::vector<std::vector<std::vector<double>>>> blocks_a(
      blocks_per_dim,
      std::vector<std::vector<std::vector<double>>>(
          blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));

  std::vector<std::vector<std::vector<std::vector<double>>>> blocks_b(
      blocks_per_dim,
      std::vector<std::vector<std::vector<double>>>(
          blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));

  std::vector<std::vector<std::vector<std::vector<double>>>> blocks_c(
      blocks_per_dim,
      std::vector<std::vector<std::vector<double>>>(
          blocks_per_dim, std::vector<std::vector<double>>(block_dim, std::vector<double>(block_dim, 0.0))));

  InitializeBlocks(source_a, source_b, blocks_a, blocks_b, block_dim, blocks_per_dim);
  RunCannonCycle(blocks_a, blocks_b, blocks_c, block_dim, blocks_per_dim);

  std::vector<std::vector<double>> result(matrix_size, std::vector<double>(matrix_size, 0.0));
  AssembleOutput(blocks_c, result, block_dim, blocks_per_dim);

  GetOutput() = std::move(result);
  return true;
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmTbb::PostProcessingImpl() {
  return true;
}

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
