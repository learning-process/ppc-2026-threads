#pragma once

#include <vector>

#include "fatehov_k_gaussian/common/include/common.hpp"
#include "task/include/task.hpp"

namespace fatehov_k_gaussian {

class FatehovKGaussianSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit FatehovKGaussianSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<float> kernel_;
  int kernel_size_ = 0;
};

}  // namespace fatehov_k_gaussian
