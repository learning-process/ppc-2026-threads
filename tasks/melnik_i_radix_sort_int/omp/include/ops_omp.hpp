#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "task/include/task.hpp"

namespace melnik_i_radix_sort_int {

class MelnikIRadixSortIntOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit MelnikIRadixSortIntOMP(const InType &in);

  struct Range {
    std::size_t begin = 0;
    std::size_t end = 0;
  };

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void CountingSortByByte(const std::vector<int> &source, std::vector<int> &destination, std::size_t begin,
                                 std::size_t end, std::int64_t exp, std::int64_t offset);
  static void RadixSortRange(std::vector<int> &data, std::vector<int> &buffer, std::size_t begin, std::size_t end);
  static void MergeRanges(const std::vector<int> &source, std::vector<int> &destination, Range left, Range right,
                          std::size_t write_begin);
  static void MergeSortedRanges(std::vector<int> &data, std::vector<int> &buffer, const std::vector<Range> &ranges);
};

}  // namespace melnik_i_radix_sort_int
