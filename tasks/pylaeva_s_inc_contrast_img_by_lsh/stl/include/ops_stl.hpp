#pragma once

#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSIncContrastImgByLshSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit PylaevaSIncContrastImgByLshSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void FindLocalMinMax(std::vector<unsigned char> &loc_mins, std::vector<unsigned char> &loc_maxs);
  void ApplyLinearStretching(int min, int max);

  int thread_count_ = 0;
  int chunk_size_ = 0;
  int size_ = 0;
};

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
