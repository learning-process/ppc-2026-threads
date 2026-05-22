#pragma once

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "task/include/task.hpp"

namespace klimovich_v_crs_complex_mat_mul {

class KlimovichVCrsComplexMatMulSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit KlimovichVCrsComplexMatMulSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static CrsMatrix MultiplyCrs(const CrsMatrix &lhs, const CrsMatrix &rhs);
};

}  // namespace klimovich_v_crs_complex_mat_mul
