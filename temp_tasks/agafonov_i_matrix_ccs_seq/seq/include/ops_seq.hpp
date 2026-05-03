#pragma once

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace agafonov_i_matrix_ccs_seq {

class AgafonovIMatrixCCSSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit AgafonovIMatrixCCSSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace agafonov_i_matrix_ccs_seq
