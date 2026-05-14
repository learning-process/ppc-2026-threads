#pragma once

#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

class ChetverikovaEShellSortSimpleMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ChetverikovaEShellSortSimpleMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &data);
  static std::vector<int> MergeSort(const std::vector<int> &left, const std::vector<int> &right);
};

}  // namespace chetverikova_e_shell_sort_simple_merge
