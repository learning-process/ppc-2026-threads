#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"

namespace terekhov_d_seq_gauss_vert {

class TerekhovDGaussVertSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit TerekhovDGaussVertSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ProcessBands(OutType &output);
  void ProcessBand(OutType &output, int padded_width, int band, int band_width);
  void ProcessPixel(OutType &output, int padded_width, int row, int col);

  int width_ = 0;
  int height_ = 0;
  static constexpr int num_bands = 4;
  std::vector<int> padded_image_;
};

}  // namespace terekhov_d_seq_gauss_vert
