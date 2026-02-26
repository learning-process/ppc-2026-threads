#pragma once

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

  int width_;
  int height_;
  std::vector<int> padded_image_;
};

}  // namespace terekhov_d_seq_gauss_vert
