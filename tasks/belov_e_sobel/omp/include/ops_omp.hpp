#pragma once

#include "belov_e_sobel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace belov_e_sobel {

class BelovESobelOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit BelovESobelOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace belov_e_sobel
