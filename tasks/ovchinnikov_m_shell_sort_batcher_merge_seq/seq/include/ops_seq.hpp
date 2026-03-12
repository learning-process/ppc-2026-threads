#pragma once

#include <vector>

#include "ovchinnikov_m_shell_sort_batcher_merge_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge_seq {

class OvchinnikovMShellSortBatcherMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit OvchinnikovMShellSortBatcherMergeSEQ(const InType &in);

 private:
  static void ShellSort(std::vector<int> &data);
  static std::vector<int> BatcherOddEvenMerge(const std::vector<int> &left, const std::vector<int> &right);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  InType data_;
};

}  // namespace ovchinnikov_m_shell_sort_batcher_merge_seq
