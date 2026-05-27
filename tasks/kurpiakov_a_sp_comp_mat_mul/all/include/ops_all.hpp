#pragma once

#include "kurpiakov_a_sp_comp_mat_mul/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kurpiakov_a_sp_comp_mat_mul {

class KurpiakovACRSMatMulALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KurpiakovACRSMatMulALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kurpiakov_a_sp_comp_mat_mul
