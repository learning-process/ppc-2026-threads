#include "safronov_m_multiplication_matrix_blockscheme_cannon/tbb/include/ops_tbb.hpp"

#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/parallel_for.h"
#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

SafronovMMultiplicationMatrixBlockSchemeCannonTBB::SafronovMMultiplicationMatrixBlockSchemeCannonTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonTBB::ValidationImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  const auto &matrix_a = std::get<1>(in);
  const auto &matrix_b = std::get<2>(in);
  return (size_block > 0) && (!matrix_a.empty() && !matrix_b.empty()) && (matrix_a.size() == matrix_a[0].size()) &&
         (matrix_b.size() == matrix_b[0].size()) && (matrix_a.size() == matrix_b.size()) &&
         (matrix_a.size() % size_block == 0);
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonTBB::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void SafronovMMultiplicationMatrixBlockSchemeCannonTBB::MultiplyingBlocks(std::vector<std::vector<double>> &block_a,
                                                                          std::vector<std::vector<double>> &block_b,
                                                                          std::vector<std::vector<double>> &block_c,
                                                                          int size_block) {
  for (int i = 0; i < size_block; i++) {
    for (int j = 0; j < size_block; j++) {
      for (int k = 0; k < size_block; k++) {
        block_c[i][j] += block_a[i][k] * block_b[k][j];
      }
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonTBB::ShiftingBlocksMatrixALeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a, int columns) {
  tbb::parallel_for(tbb::blocked_range<int>(0, columns), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); ++i) {
      std::vector<std::vector<double>> tmp = std::move(matrix_blocks_a[i][0]);
      for (int j = 1; j < columns; j++) {
        matrix_blocks_a[i][j - 1] = std::move(matrix_blocks_a[i][j]);
      }
      matrix_blocks_a[i][columns - 1] = std::move(tmp);
    }
  });
}

void SafronovMMultiplicationMatrixBlockSchemeCannonTBB::ShiftingBlocksMatrixBUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b, int columns) {
  tbb::parallel_for(tbb::blocked_range<int>(0, columns), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); ++i) {
      std::vector<std::vector<double>> tmp = std::move(matrix_blocks_b[0][i]);
      for (int j = 1; j < columns; j++) {
        matrix_blocks_b[j - 1][i] = std::move(matrix_blocks_b[j][i]);
      }
      matrix_blocks_b[columns - 1][i] = std::move(tmp);
    }
  });
}

void SafronovMMultiplicationMatrixBlockSchemeCannonTBB::AlgorithmCannon(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c, int size_block, int columns_blocks) {
  for (int i = 0; i < columns_blocks; i++) {
    tbb::parallel_for(tbb::blocked_range2d<int>(0, columns_blocks, 0, columns_blocks),
                      [&](const tbb::blocked_range2d<int> &range) {
      for (int j = range.rows().begin(); j != range.rows().end(); ++j) {
        for (int k = range.cols().begin(); k != range.cols().end(); ++k) {
          SafronovMMultiplicationMatrixBlockSchemeCannonTBB::MultiplyingBlocks(
              matrix_blocks_a[j][k], matrix_blocks_b[j][k], matrix_blocks_c[j][k], size_block);
        }
      }
    });

    if (i < columns_blocks - 1) {
      SafronovMMultiplicationMatrixBlockSchemeCannonTBB::ShiftingBlocksMatrixALeft(matrix_blocks_a, columns_blocks);
      SafronovMMultiplicationMatrixBlockSchemeCannonTBB::ShiftingBlocksMatrixBUp(matrix_blocks_b, columns_blocks);
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonTBB::FillingResultingMatrix(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c,
    std::vector<std::vector<double>> &matrix_c, int size_block, int columns_blocks) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, columns_blocks, 0, columns_blocks),
                    [&](const tbb::blocked_range2d<int> &range) {
    for (int i = range.rows().begin(); i != range.rows().end(); ++i) {
      for (int j = range.cols().begin(); j != range.cols().end(); ++j) {
        for (int k = 0; k < size_block; k++) {
          for (int col = 0; col < size_block; col++) {
            matrix_c[(i * size_block) + k][(j * size_block) + col] = matrix_blocks_c[i][j][k][col];
          }
        }
      }
    }
  });
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonTBB::RunImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  const auto &matrix_a = std::get<1>(in);
  const auto &matrix_b = std::get<2>(in);
  int n = static_cast<int>(matrix_a.size());
  int columns_blocks = n / size_block;
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_a(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_b(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_c(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block, 0.0))));

  tbb::parallel_for(tbb::blocked_range2d<int>(0, columns_blocks, 0, columns_blocks),
                    [&](const tbb::blocked_range2d<int> &range) {
    for (int i = range.rows().begin(); i != range.rows().end(); ++i) {
      for (int j = range.cols().begin(); j != range.cols().end(); ++j) {
        int shift = (i + j) % columns_blocks;
        for (int k = 0; k < size_block; k++) {
          for (int col = 0; col < size_block; col++) {
            matrix_blocks_a[i][j][k][col] = matrix_a[(i * size_block) + k][(shift * size_block) + col];
            matrix_blocks_b[i][j][k][col] = matrix_b[(shift * size_block) + k][(j * size_block) + col];
          }
        }
      }
    }
  });

  SafronovMMultiplicationMatrixBlockSchemeCannonTBB::AlgorithmCannon(matrix_blocks_a, matrix_blocks_b, matrix_blocks_c,
                                                                     size_block, columns_blocks);

  std::vector<std::vector<double>> matrix_c(n, std::vector<double>(n));
  SafronovMMultiplicationMatrixBlockSchemeCannonTBB::FillingResultingMatrix(matrix_blocks_c, matrix_c, size_block,
                                                                            columns_blocks);
  GetOutput() = std::move(matrix_c);
  return true;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonTBB::PostProcessingImpl() {
  return true;
}

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
