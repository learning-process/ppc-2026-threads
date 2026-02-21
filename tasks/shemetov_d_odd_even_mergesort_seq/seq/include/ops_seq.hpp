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
  static void CompareExchange(int &a, int &b);

  void PerfectUnshuffle(std::vector<int> &array, size_t left, size_t right);
  void PerfectShuffle(std::vector<int> &array, size_t left, size_t right);
  void Merge(std::vector<int> &array, size_t left, size_t right);
  void OddEvenMergesort(std::vector<int> &array, size_t left, size_t right);
  void SortWithPadding(std::vector<int> &array, size_t size, size_t power);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> array_;
  std::vector<int> buffer_;

  size_t size_{1};
  size_t power_{1};
};

}  // namespace shemetov_d_odd_even_mergesort
