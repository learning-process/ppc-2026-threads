#pragma once

#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

class ShkrylevaSShellMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ShkrylevaSShellMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void ShellSort(std::vector<int> &data);
  static std::vector<int> MergeSortedVectors(const std::vector<int> &left, const std::vector<int> &right);

  InType input_data_;
  OutType output_data_;
};

}  // namespace shkryleva_s_shell_sort_simple_merge
