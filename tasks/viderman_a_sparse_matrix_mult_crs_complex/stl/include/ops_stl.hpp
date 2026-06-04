#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

class VidermanASparseMatrixMultCRSComplexSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit VidermanASparseMatrixMultCRSComplexSTL(const InType &in);

 private:
  void ProcessRows(int start_row, int end_row, std::vector<std::vector<Complex>> &local_values,
                   std::vector<std::vector<int>> &local_cols) const;

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  const CRSMatrix *a_{nullptr};
  const CRSMatrix *b_{nullptr};
};

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
