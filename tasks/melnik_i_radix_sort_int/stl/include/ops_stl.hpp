#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "task/include/task.hpp"

namespace melnik_i_radix_sort_int {

constexpr int kBitsPerPass = 8;
constexpr std::size_t kBuckets = 1U << kBitsPerPass;

class MelnikIRadixSortIntSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit MelnikIRadixSortIntSTL(const InType &in);

  struct Range {
    std::size_t begin = 0;
    std::size_t end = 0;
    std::vector<int> *buffer = nullptr;
  };

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void CountingSortByByte(const std::vector<int> &source, std::vector<int> &destination, std::size_t begin,
                                 std::size_t end, std::int64_t exp, std::int64_t offset,
                                 const std::array<std::size_t, kBuckets> &count);
  static std::vector<int> *RadixSortRange(std::vector<int> &data, std::vector<int> &buffer, std::size_t begin,
                                          std::size_t end);
  static void MergeRanges(const std::vector<int> &left_source, const std::vector<int> &right_source,
                          std::vector<int> &destination, Range left, Range right, std::size_t write_begin);
  static void MergeSortedRanges(std::vector<int> &data, std::vector<int> &buffer, std::vector<Range> &ranges);

  static void EnsureRangeInBuffer(Range &range, std::vector<int> *target_buffer);
  static Range ProcessMergePair(std::size_t pair_index, const std::vector<Range> &current_ranges,
                                std::vector<int> *majority_source, std::vector<int> *level_dest);
};

}  // namespace melnik_i_radix_sort_int
