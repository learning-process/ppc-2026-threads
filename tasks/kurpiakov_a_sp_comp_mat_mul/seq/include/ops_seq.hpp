#pragma once

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

class KurpiskovACRSMatMulSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KurpiskovACRSMatMulSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kurpiakov_a_sp_comp_mat_mul
