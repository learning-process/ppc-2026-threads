#include "timur_a_cannon/seq/include/ops_seq.hpp"

#include <utility>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"

nnamespace timur_a_cannon {
  TimurACannonMatrixMultiplication::TimurACannonMatrixMultiplication(const InType &in) {
    SetTypeOfTask(GetStaticTypeOfTask());
    GetInput() = in;
  }
  bool TimurACannonMatrixMultiplication::ValidationImpl() {
    const auto &in = GetInput();
    int size_block = std::get<0>(in);
    const auto &matrix_a = std::get<1>(in);
    const auto &matrix_b = std::get<2>(in);
    return (size_block > 0) && (!matrix_a.empty() && !matrix_b.empty()) && (matrix_a.size() == matrix_a[0].size()) &&
           (matrix_b.size() == matrix_b[0].size()) && (matrix_a.size() == matrix_b.size()) &&
           (matrix_a.size() % static_cast<size_t>(size_block) == 0);
  }

  bool TimurACannonMatrixMultiplication::PreProcessingImpl() {
    GetOutput().clear();
    return true;
  }

  void TimurACannonMatrixMultiplication::MultiplyingBlocks(std::vector<std::vector<double>> & block_a,
                                                           std::vector<std::vector<double>> & block_b,
                                                           std::vector<std::vector<double>> & block_c, int size_block) {
    for (int i = 0; i < size_block; i++) {
      for (int j = 0; j < size_block; j++) {
        for (int k = 0; k < size_block; k++) {
          block_c[i][j] += block_a[i][k] * block_b[k][j];
        }
      }
    }
  }

  void TimurACannonMatrixMultiplication::ShiftingBlocksMatrixALeft(
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_a, int columns) {
    for (int i = 0; i < columns; i++) {
      std::vector<std::vector<double>> tmp = std::move(matrix_blocks_a[i][0]);
      for (int j = 1; j < columns; j++) {
        matrix_blocks_a[i][j - 1] = std::move(matrix_blocks_a[i][j]);
      }
      matrix_blocks_a[i][columns - 1] = std::move(tmp);
    }
  }

  void TimurACannonMatrixMultiplication::ShiftingBlocksMatrixBUp(
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_b, int columns) {
    for (int i = 0; i < columns; i++) {
      std::vector<std::vector<double>> tmp = std::move(matrix_blocks_b[0][i]);
      for (int j = 1; j < columns; j++) {
        matrix_blocks_b[j - 1][i] = std::move(matrix_blocks_b[j][i]);
      }
      matrix_blocks_b[columns - 1][i] = std::move(tmp);
    }
  }

  void TimurACannonMatrixMultiplication::AlgorithmCannon(
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_a,
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_b,
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_c, int size_block,
      int columns_blocks) {
    for (int i = 0; i < columns_blocks; i++) {
      for (int j = 0; j < columns_blocks; j++) {
        for (int k = 0; k < columns_blocks; k++) {
          TimurACannonMatrixMultiplication::MultiplyingBlocks(matrix_blocks_a[j][k], matrix_blocks_b[j][k],
                                                              matrix_blocks_c[j][k], size_block);
        }
      }
      if (i < columns_blocks - 1) {
        TimurACannonMatrixMultiplication::ShiftingBlocksMatrixALeft(matrix_blocks_a, columns_blocks);
        TimurACannonMatrixMultiplication::ShiftingBlocksMatrixBUp(matrix_blocks_b, columns_blocks);
      }
    }
  }

  void TimurACannonMatrixMultiplication::FillingResultingMatrix(
      std::vector<std::vector<std::vector<std::vector<double>>>> & matrix_blocks_c,
      std::vector<std::vector<double>> & matrix_c, int size_block, int columns_blocks) {
    for (int i = 0; i < columns_blocks; i++) {
      for (int j = 0; j < columns_blocks; j++) {
        for (int k = 0; k < size_block; k++) {
          for (int col = 0; col < size_block; col++) {
            matrix_c[(i * size_block) + k][(j * size_block) + col] = matrix_blocks_c[i][j][k][col];
          }
        }
      }
    }
  }

  bool TimurACannonMatrixMultiplication::RunImpl() {
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

    for (int i = 0; i < columns_blocks; i++) {
      for (int j = 0; j < columns_blocks; j++) {
        int shift = (i + j) % columns_blocks;
        for (int k = 0; k < size_block; k++) {
          for (int col = 0; col < size_block; col++) {
            matrix_blocks_a[i][j][k][col] = matrix_a[(i * size_block) + k][(shift * size_block) + col];
            matrix_blocks_b[i][j][k][col] = matrix_b[(shift * size_block) + k][(j * size_block) + col];
          }
        }
      }
    }

    TimurACannonMatrixMultiplication::AlgorithmCannon(matrix_blocks_a, matrix_blocks_b, matrix_blocks_c, size_block,
                                                      columns_blocks);

    std::vector<std::vector<double>> matrix_c(n, std::vector<double>(n));
    TimurACannonMatrixMultiplication::FillingResultingMatrix(matrix_blocks_c, matrix_c, size_block, columns_blocks);

    GetOutput() = std::move(matrix_c);
    return true;
  }

  bool TimurACannonMatrixMultiplication::PostProcessingImpl() {
    return true;
  }

}  // namespace timur_a_cannon
