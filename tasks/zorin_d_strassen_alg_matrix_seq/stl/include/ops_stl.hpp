#pragma once

#include "task/include/task.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

class ZorinDStrassenAlgMatrixSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit ZorinDStrassenAlgMatrixSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zorin_d_strassen_alg_matrix_seq
