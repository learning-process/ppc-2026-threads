#pragma once

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

  void compareExchange(std::vector<int> &a, int i, int j);
  void oddEvenMerge(std::vector<int> &a, int lo, int n, int r);
  void oddEvenMergeSort(std::vector<int> &a, int lo, int n);
};

}  // namespace khruev_a_radix_sorting_int_bather_merge
