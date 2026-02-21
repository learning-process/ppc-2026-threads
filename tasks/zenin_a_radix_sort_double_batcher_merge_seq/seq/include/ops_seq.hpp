#pragma once

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
  static void BatcherOddEvenMerge(std::vector<double> &array, int lo, int hi, int step);
  static void BatcherMergeSort(std::vector<double> &array);
};

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
