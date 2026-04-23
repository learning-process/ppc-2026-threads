#pragma once

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"
#include "task/include/task.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

class BatkovFContrastEnhLinHistStretchALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit BatkovFContrastEnhLinHistStretchALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int rank_ = 0;
  int comm_size_ = 1;
  size_t image_size_ = 0;
};

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
