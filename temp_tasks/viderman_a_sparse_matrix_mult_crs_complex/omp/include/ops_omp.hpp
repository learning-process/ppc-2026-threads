#pragma once

#include "task/include/task.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

class VidermanASparseMatrixMultCRSComplexOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit VidermanASparseMatrixMultCRSComplexOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void Multiply(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c);

  const CRSMatrix *A_{nullptr};
  const CRSMatrix *B_{nullptr};
};

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
