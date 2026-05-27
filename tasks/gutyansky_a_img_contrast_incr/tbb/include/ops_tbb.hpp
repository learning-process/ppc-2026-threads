#pragma once

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gutyansky_a_img_contrast_incr {

class GutyanskyAImgContrastIncrTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit GutyanskyAImgContrastIncrTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gutyansky_a_img_contrast_incr
