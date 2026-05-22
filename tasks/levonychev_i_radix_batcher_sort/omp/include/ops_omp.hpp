#pragma once

#include <cstddef>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit LevonychevIRadixBatcherSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void CountingSort(InType &arr, size_t byte_index);
  static void BatcherMergeIterative(std::vector<int> &arr, int start_p, int threads);
  static void BatcherCompareRange(std::vector<int> &arr, int j, int k, int p2);
};

}  // namespace levonychev_i_radix_batcher_sort
