#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "terekhov_d_gauss_vert/common/include/common.hpp"

namespace terekhov_d_gauss_vert {

class TerekhovDGaussVertSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit TerekhovDGaussVertSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int width_ = 0;
  int height_ = 0;
  std::vector<int> padded_image_;
};

}  // namespace terekhov_d_gauss_vert
