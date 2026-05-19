#pragma once

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ShkrylevaSShellMergeALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace shkryleva_s_shell_sort_simple_merge
