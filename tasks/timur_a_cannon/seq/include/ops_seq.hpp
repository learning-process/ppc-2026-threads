#pragma once
#include <vector>

#include "task/include/task.hpp"
#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

class TimurACannonMatrixMultiplication : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TimurACannonMatrixMultiplication(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void ShiftingBlocksMatrixBUp(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b,
                                      int columns_blocks);
  static void ShiftingBlocksMatrixALeft(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a,
                                        int columns_blocks);
  static void MultiplyingBlocks(std::vector<std::vector<double>> &block_a, std::vector<std::vector<double>> &block_b,
                                std::vector<std::vector<double>> &block_c, int size_block);
  static void AlgorithmCannon(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a,
                              std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b,
                              std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c,
                              int size_block, int columns_blocks);
  static void FillingResultingMatrix(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c,
                                     std::vector<std::vector<double>> &matrix_c, int size_block, int columns_blocks);
};

}  // namespace timur_a_cannon
