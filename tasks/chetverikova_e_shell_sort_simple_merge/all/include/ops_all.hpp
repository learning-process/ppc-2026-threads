#pragma once

#include <vector>

#include "chetverikova_e_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chetverikova_e_shell_sort_simple_merge {

class ChetverikovaEShellSortSimpleMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ChetverikovaEShellSortSimpleMergeALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &data);
  static std::vector<int> MergeTwoSortedVectors(const std::vector<int> &a, const std::vector<int> &b);
  static void CalculateCountsAndDisplacements(int global_size, int processes_count, std::vector<int> &counts,
                                              std::vector<int> &displacements);
  static std::vector<int> MergeLocalBuffers(std::vector<std::vector<int>> &local_buffers);
};

}  // namespace chetverikova_e_shell_sort_simple_merge
