#pragma once

#include "task/include/task.hpp"
#include "viderman_a_sparse_matrix_mult_crs_complex/common/include/common.hpp"

namespace viderman_a_sparse_matrix_mult_crs_complex {

class VidermanASparseMatrixMultCRSComplexSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit VidermanASparseMatrixMultCRSComplexSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void multiply(const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c);

  const CRSMatrix *a_;
  const CRSMatrix *b_;
};

}  // namespace viderman_a_sparse_matrix_mult_crs_complex
