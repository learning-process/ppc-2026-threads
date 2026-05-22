#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "task/include/task.hpp"

namespace melnik_i_radix_sort_int {

class MelnikIRadixSortIntALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MelnikIRadixSortIntALL(const InType &in);

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

  static void MergeSortedRangesParallel(std::vector<int> &data, std::vector<int> &buffer,
                                        const std::vector<Range> &ranges, int num_threads);

  static void MergeLevelSequential(std::vector<Range> &next, const std::vector<Range> &cur, std::vector<int> *src,
                                   std::vector<int> *dst, std::size_t pairs);

  static void MergeLevelParallel(std::vector<Range> &next, const std::vector<Range> &cur, std::vector<int> *src,
                                 std::vector<int> *dst, std::size_t pairs, int active);

  static void LocalSortChunk(std::vector<int> &local, int local_count, int local_threads, int num_threads);

  static void GlobalMergeChunks(std::vector<int> &output, int total_size, int num_ranks,
                                const std::vector<int> &send_counts, const std::vector<int> &displacements,
                                int num_threads);
};

}  // namespace melnik_i_radix_sort_int
