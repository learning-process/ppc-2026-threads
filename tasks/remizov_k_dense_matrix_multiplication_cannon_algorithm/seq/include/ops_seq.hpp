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
};

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
