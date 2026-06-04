#pragma once

#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ashihmin_d_mult_matr_crs {

class AshihminDMultMatrCrsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit AshihminDMultMatrCrsALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyRow(int global_row_idx, int local_idx, const CRSMatrix &matrix_a, const CRSMatrix &matrix_b,
                          std::vector<std::vector<int>> &local_cols, std::vector<std::vector<double>> &local_vals);
};

}  // namespace ashihmin_d_mult_matr_crs
