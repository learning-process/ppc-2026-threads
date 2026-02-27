#pragma once

#include "task/include/task.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit VasilievMShellSortBatcherMergeSEQ(const InType &in);
  static void ShellSort(std::vector<ValType> &vec);
  static std::vector<ValType> BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r);
  static std::vector<ValType> Merge(std::vector<ValType> &a, std::vector<ValType> &b);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace vasiliev_m_shell_sort_batcher_merge
