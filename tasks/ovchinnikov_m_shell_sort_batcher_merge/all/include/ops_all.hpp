#pragma once

#include <vector>

#include "ovchinnikov_m_shell_sort_batcher_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ovchinnikov_m_shell_sort_batcher_merge {

class OvchinnikovMShellSortBatcherMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit OvchinnikovMShellSortBatcherMergeALL(const InType &in);

 private:
  static void ShellSort(std::vector<int> &data);
  static std::vector<int> BatcherOddEvenMerge(const std::vector<int> &left, const std::vector<int> &right);
  static void LocalSort(std::vector<int> &local_data);
  static void TreeMerge(int rank, int active_processes, int chunk_size, std::vector<int> &local_data);

  void ScatterData(int rank, int active_processes, int chunk_size, std::size_t padded_size,
                   std::vector<int> &local_data);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace ovchinnikov_m_shell_sort_batcher_merge
