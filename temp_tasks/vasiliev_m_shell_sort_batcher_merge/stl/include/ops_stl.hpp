#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

class VasilievMShellSortBatcherMergeSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit VasilievMShellSortBatcherMergeSTL(const InType &in);
  static std::vector<size_t> ChunkBoundaries(size_t vec_size, int threads);
  static void ShellSort(std::vector<ValType> &vec, std::vector<size_t> &bounds, size_t num_threads);
  static void CycleMerge(std::vector<ValType> &vec, std::vector<ValType> &buffer, std::vector<size_t> &bounds,
                         size_t size, size_t num_threads);
  static std::vector<ValType> BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r);
  static void SplitEvenOdd(std::vector<ValType> &vec, std::vector<ValType> &even, std::vector<ValType> &odd);
  static std::vector<ValType> Merge(std::vector<ValType> &a, std::vector<ValType> &b);
  static void ShellSortChunk(std::vector<ValType> &vec, size_t first, size_t last);
  static void MergeChunks(std::vector<ValType> &vec, std::vector<ValType> &buffer, size_t start, size_t middle,
                          size_t end_pos);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace vasiliev_m_shell_sort_batcher_merge
