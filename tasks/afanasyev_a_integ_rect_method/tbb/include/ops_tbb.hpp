#pragma once

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace afanasyev_a_integ_rect_method {

class AfanasyevAIntegRectMethodTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit AfanasyevAIntegRectMethodTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace afanasyev_a_integ_rect_method
