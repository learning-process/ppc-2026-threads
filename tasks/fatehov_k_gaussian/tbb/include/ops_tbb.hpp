#pragma once

#include "fatehov_k_gaussian/common/include/common.hpp"
#include "task/include/task.hpp"

namespace fatehov_k_gaussian {

class FatehovKGaussianTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit FatehovKGaussianTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace fatehov_k_gaussian
