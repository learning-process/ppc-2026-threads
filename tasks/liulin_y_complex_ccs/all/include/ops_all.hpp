#pragma once

#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace liulin_y_complex_ccs {

class LiulinYComplexCcsAll : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LiulinYComplexCcsAll(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace liulin_y_complex_ccs
