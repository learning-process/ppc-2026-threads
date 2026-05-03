#pragma once

#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

class ChetverikovaEShellSortSimpleMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ChetverikovaEShellSortSimpleMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &data);
};

}  // namespace chetverikova_e_shell_sort_simple_merge
