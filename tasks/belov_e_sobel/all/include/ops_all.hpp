#pragma once

#include "belov_e_sobel/common/include/common.hpp"
#include "task/include/task.hpp"

namespace belov_e_sobel {

class BelovESobelALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit BelovESobelALL(const InType &in);

 private:
  bool ValidationImpl() override;
  static bool PreProcessingImpl() override;
  static bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace belov_e_sobel
