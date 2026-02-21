#pragma once

#include "kapanova_s_increase_contrast_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kapanova_s_increase_contrast_seq {

class KapanovaSIncreaseContrastSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KapanovaSIncreaseContrastSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  uint8_t min_pixel_;
  uint8_t max_pixel_;
};

}  // namespace kapanova_s_increase_contrast_seq
