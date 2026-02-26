#pragma once

#include <cstddef>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shemetov_d_odd_even_mergesort {

class ShemetovDOddEvenMergeSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit ShemetovDOddEvenMergeSortSEQ(const InType &in);

 private:
  void RadixSort(std::vector<int> &array, size_t left, size_t right);
  static void OddEvenMerge(std::vector<int> &array, size_t start, size_t segment);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> array_;
  std::vector<int> buffer_;
  std::vector<int> position_;

  int offset_{0};
  size_t size_{0};
  size_t power_{0};
};

}  // namespace shemetov_d_odd_even_mergesort
