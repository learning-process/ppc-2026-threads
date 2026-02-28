#pragma once

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace leonova_a_radix_merge_sort {

class LeonovaARadixMergeSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit LeonovaARadixMergeSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixMergeSort(std::vector<int64_t> &arr, size_t left, size_t right);
  static void RadixSort(std::vector<int64_t> &arr, size_t left, size_t right);
  static void SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right);

  static uint64_t ToUnsigned(int64_t value);

  static constexpr size_t kRadixThreshold = 32;  // Порог для переключения на RadixSort
  static constexpr int kByteSize = 8;
  static constexpr int kNumBytes = 8;
  static constexpr int kNumCounters = 256;
};

}  // namespace leonova_a_radix_merge_sort
