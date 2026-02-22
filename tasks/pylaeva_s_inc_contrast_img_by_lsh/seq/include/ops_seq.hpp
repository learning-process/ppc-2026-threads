#pragma once

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSIncContrastImgByLshSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PylaevaSIncContrastImgByLshSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
