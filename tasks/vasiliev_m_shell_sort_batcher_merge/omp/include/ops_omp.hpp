#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergeOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit VasilievMShellSortBatcherMergeOMP(const InType &in);
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
