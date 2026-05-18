#pragma once

#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ShkrylevaSShellMergeALL(const InType &in);

  static std::vector<int> SimpleMerge(const std::vector<int> &a, const std::vector<int> &b);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &arr, int left, int right);
  static void Merge(std::vector<int> &arr, int left, int mid, int right, std::vector<int> &buffer);
  static void ParallelShellSort(std::vector<int> &arr);
};

}  // namespace shkryleva_s_shell_sort_simple_merge
