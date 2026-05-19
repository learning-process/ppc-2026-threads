#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "task/include/task.hpp"
#include "tbb/tbb.h"

namespace leonova_a_radix_merge_sort {

class LeonovaARadixMergeSortTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit LeonovaARadixMergeSortTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixMergeSort(std::vector<int64_t> &arr, size_t left, size_t right);
  static void RadixSort(std::vector<int64_t> &arr, size_t left, size_t right);
  static void SequentialRadixSort(std::vector<int64_t> &arr, size_t left, size_t right);
  static void SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right);

  static uint64_t ToUnsignedValue(int64_t value);

  static std::pair<size_t, size_t> GetChunk(size_t tid, size_t num_threads, size_t size);

  using CounterRow = std::vector<size_t>;

  static void ScatterParallel(const std::vector<uint64_t> &keys, const std::vector<int64_t> &arr, size_t left,
                              size_t size, int shift, std::vector<CounterRow> &local_offsets,
                              std::vector<int64_t> &temp_arr, std::vector<uint64_t> &temp_keys, size_t num_threads);

  static void BuildOffsets(const std::vector<CounterRow> &local_counts, std::vector<CounterRow> &local_offsets,
                           CounterRow &global_counts, size_t num_threads);

  static void ReduceCounts(const std::vector<CounterRow> &local_counts, CounterRow &global_counts);

  static void CountBytesParallel(const std::vector<uint64_t> &keys, size_t size, int shift,
                                 std::vector<CounterRow> &local_counts, size_t num_threads);

  static void FillUnsignedKeys(const std::vector<int64_t> &arr, size_t left, size_t size, std::vector<uint64_t> &keys,
                               size_t num_threads);

  static constexpr size_t kRadixThreshold = 131072;
  static constexpr size_t kMinParallelSize = 10000;
  static constexpr int kByteSize = 8;
  static constexpr int kNumBytes = 8;
  static constexpr int kNumCounters = 256;
  static constexpr uint64_t kSignBitMask = 0x8000000000000000ULL;

  static tbb::task_arena &GetTbbArena();
};

}  // namespace leonova_a_radix_merge_sort
