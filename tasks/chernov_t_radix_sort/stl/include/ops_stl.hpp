#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace chernov_t_radix_sort {

class ChernovTRadixSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ChernovTRadixSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSortLSDSequential(std::vector<int> &data);
  static void RadixSortLSDParallel(std::vector<int> &data, int num_threads);

  static void ConvertToIntegers(const std::vector<int> &input, std::vector<uint32_t> &output, size_t start, size_t end,
                                uint32_t sign_mask);
  static void ConvertFromIntegers(const std::vector<uint32_t> &input, std::vector<int> &output, size_t start,
                                  size_t end, uint32_t sign_mask);
  static void ComputeLocalHistograms(const std::vector<uint32_t> &data, std::vector<std::vector<int>> &local_counts,
                                     size_t start, size_t end, int shift, int thread_idx);
  static void ComputeGlobalStarts(const std::vector<std::vector<int>> &local_counts, std::vector<int> &global_start,
                                  int k_radix, int num_threads);
  static void ComputeThreadOffsets(const std::vector<std::vector<int>> &local_counts,
                                   std::vector<std::vector<int>> &thread_offset, int k_radix, int num_threads);
  static void ScatterElements(const std::vector<uint32_t> &input, std::vector<uint32_t> &output,
                              const std::vector<int> &global_start, const std::vector<std::vector<int>> &thread_offset,
                              std::vector<std::vector<int>> &local_counter, size_t start, size_t end, int shift,
                              int thread_idx);
};

}  // namespace chernov_t_radix_sort
