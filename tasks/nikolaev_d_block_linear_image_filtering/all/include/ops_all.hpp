#pragma once

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"
#include "task/include/task.hpp"

namespace nikolaev_d_block_linear_image_filtering {

class NikolaevDBlockLinearImageFilteringALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit NikolaevDBlockLinearImageFilteringALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace nikolaev_d_block_linear_image_filtering
