#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "task/include/task.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

class ZeninARadixSortDoubleBatcherMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ZeninARadixSortDoubleBatcherMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void BlocksComparing(std::vector<double> &arr, size_t i, size_t j);
  static uint64_t PackDouble(double v) noexcept;
  static double UnpackDouble(uint64_t k) noexcept;
  static void LSDRadixSort(std::vector<double> &arr);
  static void BatcherOddEvenMerge(std::vector<double> &arr, size_t n);
};

}  // namespace zenin_a_radix_sort_double_batcher_merge
