#pragma once

#include <vector>

#include "sakharov_a_shell_sorting_with_merging_butcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sakharov_a_shell_sorting_with_merging_butcher {

class SakharovAShellButcherSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SakharovAShellButcherSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &data);
  static std::vector<int> MergeSortedVectors(const std::vector<int> &left, const std::vector<int> &right);
  static std::vector<int> BatcherOddEvenMerge(const std::vector<int> &left, const std::vector<int> &right);

  InType input_data_;
  OutType output_data_;
};

}  // namespace sakharov_a_shell_sorting_with_merging_butcher
