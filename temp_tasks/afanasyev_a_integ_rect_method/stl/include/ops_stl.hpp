#pragma once

#include "afanasyev_a_integ_rect_method/common/include/common.hpp"
#include "task/include/task.hpp"

namespace afanasyev_a_integ_rect_method {

class AfanasyevAIntegRectMethodSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit AfanasyevAIntegRectMethodSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace afanasyev_a_integ_rect_method
