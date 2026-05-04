#pragma once

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dergynov_s_integrals_multistep_rectangle {

class DergynovSIntegralsMultistepRectangleALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit DergynovSIntegralsMultistepRectangleALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace dergynov_s_integrals_multistep_rectangle
