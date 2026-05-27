#pragma once

#include <cstddef>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"
#include "task/include/task.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

class PylaevaSIncContrastImgByLshALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit PylaevaSIncContrastImgByLshALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType local_data_;
  OutType local_out_;
  size_t local_size_{};
  std::vector<int> recv_displs_;
  std::vector<int> recv_counts_;
};

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
