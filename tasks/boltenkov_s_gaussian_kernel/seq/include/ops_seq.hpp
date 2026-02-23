#pragma once

#include <vector>

#include "boltenkov_s_gaussian_kernel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace boltenkov_s_gaussian_kernel {

class BoltenkovSGaussianKernelSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit BoltenkovSGaussianKernelSEQ(const InType &in);

 private:
  std::vector<std::vector<int>> kernel_;
  int shift_ = 4;
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace boltenkov_s_gaussian_kernel
