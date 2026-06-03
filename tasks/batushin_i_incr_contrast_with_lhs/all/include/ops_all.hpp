#pragma once

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace batushin_i_incr_contrast_with_lhs {

class BatushinIIncrContrastWithLhsALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }

  explicit BatushinIIncrContrastWithLhsALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace batushin_i_incr_contrast_with_lhs
