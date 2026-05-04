#pragma once

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

class GoriachevaKMultSparseComplexMatrixCcsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit GoriachevaKMultSparseComplexMatrixCcsALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
