#pragma once

#include "dorogin_v_bin_img_conv_hull_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull_all {

class DoroginVConvHullAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit DoroginVConvHullAll(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType local_out_;
};

}  // namespace dorogin_v_bin_img_conv_hull_all
