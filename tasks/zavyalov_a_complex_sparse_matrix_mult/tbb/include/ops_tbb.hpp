#pragma once

#include "task/include/task.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {  // comm for retrigerring CI

class ZavyalovAComplSparseMatrMultTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ZavyalovAComplSparseMatrMultTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static SparseMatrix MultiplicateWithTbb(const SparseMatrix &matr_a, const SparseMatrix &matr_b);
};

}  // namespace zavyalov_a_compl_sparse_matr_mult
