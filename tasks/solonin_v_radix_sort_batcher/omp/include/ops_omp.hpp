#pragma once
#include <cstddef>
#include <vector>

#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "task/include/task.hpp"

namespace solonin_v_radix_sort_batcher {
class RadixSortBatcherOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit RadixSortBatcherOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SortByDigit(std::vector<int> &data, size_t pos);
  static void SortChunk(std::vector<int> &data, int left, int right);
  static void BatcherNetwork(std::vector<int> &data, int p, int nthreads);
  static void CompareSwapRange(std::vector<int> &data, int offset, int step, int half);
};
}  // namespace solonin_v_radix_sort_batcher
