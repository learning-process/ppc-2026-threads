#pragma once

#include "kazennova_a_sobel_operator/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kazennova_a_sobel_operator {

class SobelSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SobelSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kazennova_a_sobel_operator
