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
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<int> array_;

  size_t size_{1};
  size_t power_{1};
};

}  // namespace shemetov_d_odd_even_mergesort
