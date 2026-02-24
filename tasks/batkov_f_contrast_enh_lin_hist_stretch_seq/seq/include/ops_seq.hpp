#pragma once

#include "batkov_f_contrast_enh_lin_hist_stretch_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch_seq {

class BatkovFContrastEnhLinHistStretchSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BatkovFContrastEnhLinHistStretchSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace batkov_f_contrast_enh_lin_hist_stretch_seq
