#pragma once

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SakharovAShellButcherOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
