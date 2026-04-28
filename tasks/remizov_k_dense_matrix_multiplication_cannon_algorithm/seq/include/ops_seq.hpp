#pragma once

#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

class RemizovKDenseMatrixMultiplicationCannonAlgorithm : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit RemizovKDenseMatrixMultiplicationCannonAlgorithm(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyBlock(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b,
                            std::vector<std::vector<double>> &c, int block_size);

  static void ShiftBlocksLeft(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks,
                              int block_count);

  static void ShiftBlocksUp(std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks, int block_count);

  static void RunCannonCycle(std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
                             std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks,
                             std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks, int block_size,
                             int block_count);

  static void AssembleOutput(std::vector<std::vector<std::vector<std::vector<double>>>> &c_blocks,
                             std::vector<std::vector<double>> &output, int block_size, int block_count);

  static void InitializeBlocks(const std::vector<std::vector<double>> &matrix_a,
                               const std::vector<std::vector<double>> &matrix_b,
                               std::vector<std::vector<std::vector<std::vector<double>>>> &a_blocks,
                               std::vector<std::vector<std::vector<std::vector<double>>>> &b_blocks, int block_size,
                               int block_count);
};

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
