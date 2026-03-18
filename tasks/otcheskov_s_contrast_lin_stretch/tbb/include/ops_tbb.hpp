#pragma once

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace otcheskov_s_contrast_lin_stretch {

class OtcheskovSContrastLinStretchTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit OtcheskovSContrastLinStretchTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace otcheskov_s_contrast_lin_stretch
