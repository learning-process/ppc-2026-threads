#pragma once

#include "fedoseev_linear_image_filtering_vertical/common/include/common.hpp"
#include "task/include/task.hpp"

namespace fedoseev_linear_image_filtering_vertical {

class LinearImageFilteringVerticalAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LinearImageFilteringVerticalAll(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType local_out_;
};

}  // namespace fedoseev_linear_image_filtering_vertical
