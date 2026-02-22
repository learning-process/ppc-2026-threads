#include "safronov_m_multiplication_matrix_blockscheme_cannon/seq/include/ops_seq.hpp"

#include <vector>

#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"

namespace safronov_m_multiplication_matrix_blockscheme_cannon {

SafronovMMultiplicationMatrixBlockSchemeCannon::SafronovMMultiplicationMatrixBlockSchemeCannon(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannon::ValidationImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  const auto &matrixA = std::get<1>(in);
  const auto &matrixB = std::get<2>(in);
  return (size_block > 0) && (!matrixA.empty() && !matrixB.empty()) && (matrixA.size() == matrixA[0].size()) &&
         (matrixB.size() == matrixB[0].size()) && (matrixA.size() == matrixB.size()) &&
         (matrixA.size() % size_block == 0);
}

bool SafronovMMultiplicationMatrixBlockSchemeCannon::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void SafronovMMultiplicationMatrixBlockSchemeCannon::MultiplyingBlocks(std::vector<std::vector<double>> &blockA,
                                                                       std::vector<std::vector<double>> &blockB,
                                                                       std::vector<std::vector<double>> &blockC,
                                                                       int size_block) {
  for (int i = 0; i < size_block; i++) {
    for (int j = 0; j < size_block; j++) {
      for (int k = 0; k < size_block; k++) {
        blockC[i][j] += blockA[i][k] * blockB[k][j];
      }
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannon::ShiftingBlocksMatrixALeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksA, int columns) {
  for (int i = 0; i < columns; i++) {
    std::vector<std::vector<double>> tmp = std::move(matrix_blocksA[i][0]);
    for (int j = 1; j < columns; j++) {
      matrix_blocksA[i][j - 1] = std::move(matrix_blocksA[i][j]);
    }
    matrix_blocksA[i][columns - 1] = std::move(tmp);
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannon::ShiftingBlocksMatrixBUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksB, int columns) {
  for (int i = 0; i < columns; i++) {
    std::vector<std::vector<double>> tmp = std::move(matrix_blocksB[0][i]);
    for (int j = 1; j < columns; j++) {
      matrix_blocksB[j - 1][i] = std::move(matrix_blocksB[j][i]);
    }
    matrix_blocksB[columns - 1][i] = std::move(tmp);
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannon::AlgorithmCannon(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksA,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksB,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksC, int size_block, int columns_blocks) {
  for (int i = 0; i < columns_blocks; i++) {
    for (int j = 0; j < columns_blocks; j++) {
      for (int k = 0; k < columns_blocks; k++) {
        MultiplyingBlocks(matrix_blocksA[j][k], matrix_blocksB[j][k], matrix_blocksC[j][k], size_block);
      }
    }
    if (i < columns_blocks - 1) {
      ShiftingBlocksMatrixALeft(matrix_blocksA, columns_blocks);
      ShiftingBlocksMatrixBUp(matrix_blocksB, columns_blocks);
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannon::FillingResultingMatrix(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksC,
    std::vector<std::vector<double>> &matrixC, int size_block, int columns_blocks) {
  for (int i = 0; i < columns_blocks; i++) {
    for (int j = 0; j < columns_blocks; j++) {
      for (int k = 0; k < size_block; k++) {
        for (int l = 0; l < size_block; l++) {
          matrixC[i * size_block + k][j * size_block + l] = matrix_blocksC[i][j][k][l];
        }
      }
    }
  }
}

bool SafronovMMultiplicationMatrixBlockSchemeCannon::RunImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  auto &matrixA = std::get<1>(in);
  auto &matrixB = std::get<2>(in);
  int N = matrixA.size();
  int columns_blocks = N / size_block;
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocksA(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocksB(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocksC(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block, 0.0))));

  for (int i = 0; i < columns_blocks; i++) {
    for (int j = 0; j < columns_blocks; j++) {
      int shift = (i + j) % columns_blocks;
      for (int k = 0; k < size_block; k++) {
        for (int l = 0; l < size_block; l++) {
          matrix_blocksA[i][j][k][l] = matrixA[i * size_block + k][shift * size_block + l];
          matrix_blocksB[i][j][k][l] = matrixB[shift * size_block + k][j * size_block + l];
        }
      }
    }
  }
  AlgorithmCannon(matrix_blocksA, matrix_blocksB, matrix_blocksC, size_block, columns_blocks);

  std::vector<std::vector<double>> matrixC(N, std::vector<double>(N));
  FillingResultingMatrix(matrix_blocksC, matrixC, size_block, columns_blocks);
  GetOutput() = std::move(matrixC);
  return true;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannon::PostProcessingImpl() {
  return true;
}

}  // namespace safronov_m_multiplication_matrix_blockscheme_cannon
