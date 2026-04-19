#pragma once
#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace luzan_e_double_sparse_matrix_mult {

class LuzanEDoubleSparseMatrixMultOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit LuzanEDoubleSparseMatrixMultOMP(const InType &in);

  SparseMatrix CalcProdOMP(const SparseMatrix &a, const SparseMatrix &b);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace luzan_e_double_sparse_matrix_mult
