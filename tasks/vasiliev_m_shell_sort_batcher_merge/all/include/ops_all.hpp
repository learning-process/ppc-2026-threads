#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit VasilievMShellSortBatcherMergeALL(const InType &in);
  static void CalcCountsAndDispls(int n, int process_count, std::vector<int> &counts, std::vector<int> &displs);
  static std::vector<size_t> BoundsFromCounts(const std::vector<int> &counts);
  static std::vector<size_t> ChunkBoundaries(size_t vec_size, int threads);
  static void ShellSort(std::vector<ValType> &vec, std::vector<size_t> &bounds);
  static void CycleMerge(std::vector<ValType> &vec, std::vector<ValType> &buffer, std::vector<size_t> &bounds,
                         size_t size);
  static std::vector<ValType> BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r);
  static void SplitEvenOdd(std::vector<ValType> &vec, std::vector<ValType> &even, std::vector<ValType> &odd);
  static std::vector<ValType> Merge(std::vector<ValType> &a, std::vector<ValType> &b);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace vasiliev_m_shell_sort_batcher_merge
