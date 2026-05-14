#pragma once

#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit LevonychevIRadixBatcherSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void RadixSortSequential(std::vector<int> &arr);
  static void MergeAndSplit(std::vector<int> &left_block, std::vector<int> &right_block);
  static int GetNumBlocks(int n);
  static std::vector<std::vector<int>> DistributeData(const std::vector<int> &data, int num_blocks);
  static void ParallelRadixPhase(std::vector<std::vector<int>> &blocks);
  static void BatcherMergePhase(std::vector<std::vector<int>> &blocks);
  static void BatcherMergeStep(std::vector<std::vector<int>> &blocks, int p, int k);
};

}  // namespace levonychev_i_radix_batcher_sort
