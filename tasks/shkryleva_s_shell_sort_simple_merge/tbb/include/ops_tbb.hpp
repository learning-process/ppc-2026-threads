#pragma once

#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeTBB : public BaseTask {
 public:
  explicit ShkrylevaSShellMergeTBB(const InType &in);

  bool PreProcessingImpl() override;
  bool ValidationImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

 private:
  static void ShellSort(int left, int right, std::vector<int> &arr);
  static void Merge(int left, int mid, int right, std::vector<int> &arr, std::vector<int> &buffer);

  InType input_data_;
  OutType output_data_;
};

}  // namespace shkryleva_s_shell_sort_simple_merge
