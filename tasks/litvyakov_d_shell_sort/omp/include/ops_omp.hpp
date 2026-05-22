#pragma once

#include <cstdio>
#include <vector>

#include "litvyakov_d_shell_sort/common/include/common.hpp"
#include "task/include/task.hpp"
namespace litvyakov_d_shell_sort {

class LitvyakovDShellSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit LitvyakovDShellSortOMP(const InType &in);

  static void BaseShellSort(std::vector<int>::iterator first, std::vector<int>::iterator last);
  static void ShellSortMerge(std::vector<int> &vec);
  static std::vector<std::size_t> GetBounds(std::size_t n, std::size_t parts);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace litvyakov_d_shell_sort
