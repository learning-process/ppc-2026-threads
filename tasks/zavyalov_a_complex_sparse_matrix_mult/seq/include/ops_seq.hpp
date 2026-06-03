#pragma once

#include "task/include/task.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZavyalovAComplSparseMatrMultSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZavyalovAComplSparseMatrMultSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zavyalov_a_compl_sparse_matr_mult
