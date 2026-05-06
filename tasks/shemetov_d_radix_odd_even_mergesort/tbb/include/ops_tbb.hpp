#pragma once

#include <cstddef>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

class ShemetovDRadixOddEvenMergeSortTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit ShemetovDRadixOddEvenMergeSortTBB(const InType &in);

 private:
  static void RadixSort(std::vector<int> &array, size_t left, size_t right);
  static void OddEvenMerge(std::vector<int> &array, size_t start, size_t segment);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> array_;

  int offset_{0};
  size_t size_{0};
  size_t power_{0};
};

}  // namespace shemetov_d_radix_odd_even_mergesort
