#pragma once

#include "nazyrov_a_multidim_integral_rectangle/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nazyrov_a_multidim_integral_rectangle {

class NazyrovAMultidimIntegralRectangleTbb : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit NazyrovAMultidimIntegralRectangleTbb(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nazyrov_a_multidim_integral_rectangle
