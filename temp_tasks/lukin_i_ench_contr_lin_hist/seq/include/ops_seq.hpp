#pragma once

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lukin_i_ench_contr_lin_hist {

class LukinITestTaskSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LukinITestTaskSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace lukin_i_ench_contr_lin_hist
