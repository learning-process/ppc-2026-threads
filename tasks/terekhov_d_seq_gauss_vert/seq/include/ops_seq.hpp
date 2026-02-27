#pragma once

#include <vector>

#include "terekhov_d_seq_gauss_vert/common/include/common.hpp"
#include "task/include/task.hpp"

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
  
  int width_ = 0;
  int height_ = 0;
  std::vector<int> padded_image_;
};

}  // namespace terekhov_d_seq_gauss_vert
