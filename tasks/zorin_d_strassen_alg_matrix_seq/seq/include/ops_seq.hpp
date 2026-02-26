#pragma once

#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

class ZorinDStrassenAlgMatrixSEQ  : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZorinDStrassenAlgMatrixSEQ (const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace zorin_d_strassen_alg_matrix_seq
