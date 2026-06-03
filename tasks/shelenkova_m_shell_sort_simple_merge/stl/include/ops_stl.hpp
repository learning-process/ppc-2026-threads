#pragma once

#include "shelenkova_m_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shelenkova_m_shell_sort_simple_merge {

class ShelenkovaMShellSortSimpleMergeSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ShelenkovaMShellSortSimpleMergeSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace shelenkova_m_shell_sort_simple_merge
