#pragma once

#include <cstddef>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

class KhruevARadixSortingIntBatherMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KhruevARadixSortingIntBatherMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void CompareExchange(std::vector<int> &a, size_t i, size_t j);
  static void RadixSort(std::vector<int> &arr);
  static void OddEvenMerge(std::vector<int> &a, size_t n);
};

}  // namespace khruev_a_radix_sorting_int_bather_merge
