#pragma once

#include <vector>

#include "litvyakov_d_shell_sort/common/include/common.hpp"
#include "task/include/task.hpp"
namespace litvyakov_d_shell_sort {

class LitvyakovDShellSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LitvyakovDShellSortSEQ(const InType &in);

  static void BaseShellSort(std::vector<int> &vec);
  static void ShellSortMerge(std::vector<int> &vec);
  static void Merge(std::vector<int> &left, const std::vector<int> &right, std::vector<int> &vec);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace litvyakov_d_shell_sort
