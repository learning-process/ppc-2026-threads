#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

class TsyplakovKTestTaskSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit TsyplakovKTestTaskSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static std::vector<double> MultiplyRowByMatrix(const std::vector<double> &row_values,
                                                 const std::vector<int> &row_cols, const SparseMatrixCRS &b,
                                                 int &result_nnz);
};

}  // namespace tsyplakov_k_mul_double_crs_matrix
