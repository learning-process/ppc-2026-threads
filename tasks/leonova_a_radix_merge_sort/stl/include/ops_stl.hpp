#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace leonova_a_radix_merge_sort {

class LeonovaARadixMergeSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit LeonovaARadixMergeSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static uint64_t ToUnsignedValue(int64_t value);

  static void FillKeys(const std::vector<int64_t> &arr, size_t left, size_t size, std::vector<uint64_t> &keys);

  static void CountBytes(const std::vector<uint64_t> &keys, size_t begin, size_t end, int shift,
                         std::vector<size_t> &counts);

  static void ScatterBytes(const std::vector<uint64_t> &src_keys, const std::vector<int64_t> &src_arr, size_t left,
                           size_t begin, size_t end, int shift, std::vector<size_t> &offsets,
                           std::vector<int64_t> &dst_arr, std::vector<uint64_t> &dst_keys);

  static void RadixSort(std::vector<int64_t> &arr, size_t left, size_t right);

  static void SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right);

  static void RadixMergeSort(std::vector<int64_t> &arr, size_t left, size_t right);

  static constexpr int kByteSize = 8;
  static constexpr int kNumBytes = 8;
  static constexpr int kNumCounters = 256;
  static constexpr uint64_t kSignBitMask = 0x8000000000000000ULL;
};

}  // namespace leonova_a_radix_merge_sort
