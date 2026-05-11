#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

class ZeninARadixSortDoubleBatcherMergeALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ZeninARadixSortDoubleBatcherMergeALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static uint64_t PackDouble(double v) noexcept;
  static double UnpackDouble(uint64_t k) noexcept;
  static void LSDRadixSort(std::vector<double> &array);
  static void BlocksComparing(std::vector<double> &arr, size_t i, size_t step);
  static void BatcherOddEvenMerge(std::vector<double> &arr, size_t n);
  static void SortChunk(std::vector<double> &chunk, size_t chunk_size);
  void FinalMerge(size_t chunk_size, size_t pow2, size_t original_size);

  std::vector<double> local_data_;
};

}  // namespace zenin_a_radix_sort_double_batcher_merge
