#pragma once

#include "kapanova_s_sparse_matrix_mult_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs {

class KapanovaSSparseMatrixMultCCSTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit KapanovaSSparseMatrixMultCCSTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kapanova_s_sparse_matrix_mult_ccs
