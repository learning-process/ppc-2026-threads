#pragma once

#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ashihmin_d_mult_matr_crs {

class AshihminDMultMatrCrsSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit AshihminDMultMatrCrsSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyRow(int row_idx, const CRSMatrix &matrix_a, const CRSMatrix &matrix_b, std::vector<int> &row_cols,
                          std::vector<double> &row_vals);
};

}  // namespace ashihmin_d_mult_matr_crs
