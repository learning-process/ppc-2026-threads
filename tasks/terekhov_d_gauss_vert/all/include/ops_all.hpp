#pragma once

#include "task/include/task.hpp"
#include "terekhov_d_gauss_vert/common/include/common.hpp"

namespace terekhov_d_gauss_vert {

class TerekhovDGaussVertALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit TerekhovDGaussVertALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType local_out_;
};

}  // namespace terekhov_d_gauss_vert
