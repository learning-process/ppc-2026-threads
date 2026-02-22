#pragma once

#include <vector>

#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"
#include "task/include/task.hpp"

namespace safronov_m_multiplication_matrix_blockscheme_cannon {

class SafronovMMultiplicationMatrixBlockSchemeCannon : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SafronovMMultiplicationMatrixBlockSchemeCannon(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  void ShiftingBlocksMatrixBUp(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksB,
                               int columns_blocks);
  void ShiftingBlocksMatrixALeft(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksA,
                                 int columns_blocks);
  void MultiplyingBlocks(std::vector<std::vector<double>> &blockA, std::vector<std::vector<double>> &blockB,
                         std::vector<std::vector<double>> &blockC, int size_block);
  void AlgorithmCannon(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksA,
                       std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksB,
                       std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksC, int size_block,
                       int columns_blocks);
  void FillingResultingMatrix(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocksC,
                              std::vector<std::vector<double>> &matrixC, int size_block, int columns_blocks);
};

}  // namespace safronov_m_multiplication_matrix_blockscheme_cannon
