#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

class ZeninARadixSortDoubleBatcherMergeSeqseq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ZeninARadixSortDoubleBatcherMergeSeqseq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void LSDRadixSort(std::vector<double> &array);
  static uint64_t PackDouble(double v) noexcept;
  static double UnpackDouble(uint64_t k) noexcept;
  static void BatcherOddEvenMerge(std::vector<double> &array, size_t n);
  static void BatcherMergeSort(std::vector<double> &array);
  static void BlocksComparing(std::vector<double> &arr, size_t i, size_t step);
};

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
