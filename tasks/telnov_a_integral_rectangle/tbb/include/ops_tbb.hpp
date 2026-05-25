#pragma once

#include "../../modules/task/include/task.hpp"
#include "telnov_a_integral_rectangle/common/include/common.hpp"

namespace telnov_a_integral_rectangle {

class TelnovAIntegralRectangleTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit TelnovAIntegralRectangleTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace telnov_a_integral_rectangle
