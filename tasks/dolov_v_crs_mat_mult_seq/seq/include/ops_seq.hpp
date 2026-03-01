#pragma once

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dolov_v_crs_mat_mult_seq {

class DolovVCrsMatMultSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit DolovVCrsMatMultSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  [[nodiscard]] SparseMatrix TransposeMatrix(const SparseMatrix &matrix) const;
};

}  // namespace dolov_v_crs_mat_mult_seq
