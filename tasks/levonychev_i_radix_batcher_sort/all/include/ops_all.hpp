#pragma once

#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LevonychevIRadixBatcherSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void LocalRadixSort(std::vector<int> &arr);
  static void LocalParallelMerge(std::vector<std::vector<int>> &blocks);

  static void NetworkMergeAndSplit(std::vector<int> &local_data, int partner, bool keep_low);

  static void CalculateDistribution(int total_n, int size, std::vector<int> &counts, std::vector<int> &displs);
  static void LocalSortPhase(std::vector<int> &local_data);
  static void LocalBatcherMerge(std::vector<std::vector<int>> &blocks);
  static void BatcherStep(std::vector<std::vector<int>> &blocks, int pr, int k);
  static void CompareAndMergeBlocks(std::vector<int> &b1, std::vector<int> &b2);
  static void GlobalCompareExchange(std::vector<int> &local_data, int rank, int i1, int i2);
  static void GlobalBatcherStep(std::vector<int> &local_data, int rank, int size, int p, int k);
  static void GlobalSortPhase(std::vector<int> &local_data, int rank, int size);
};

}  // namespace levonychev_i_radix_batcher_sort
