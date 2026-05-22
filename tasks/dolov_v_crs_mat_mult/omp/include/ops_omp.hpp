#pragma once

#include "dolov_v_crs_mat_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dolov_v_crs_mat_mult {

class DolovVCrsMatMultOmp : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit DolovVCrsMatMultOmp(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] static SparseMatrix TransposeMatrix(const SparseMatrix &matrix);
  static double DotProduct(const SparseMatrix &matrix_a, int row_a, const SparseMatrix &matrix_b_t, int row_b);
};

}  // namespace dolov_v_crs_mat_mult
