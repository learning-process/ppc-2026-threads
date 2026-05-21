#pragma once

#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "task/include/task.hpp"

namespace klimovich_v_crs_complex_mat_mul {

class KlimovichVCrsComplexMatMulAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit KlimovichVCrsComplexMatMulAll(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void BroadcastOperand(CrsMatrix &m, int root);
  static void ComputeLocalRows(const CrsMatrix &lhs, const CrsMatrix &rhs, int row_begin, int row_end,
                               std::vector<int> &local_nnz_per_row, std::vector<int> &local_cols,
                               std::vector<Cplx> &local_vals);
};

}  // namespace klimovich_v_crs_complex_mat_mul
