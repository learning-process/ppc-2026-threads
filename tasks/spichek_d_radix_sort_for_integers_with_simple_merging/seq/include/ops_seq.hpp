#pragma once

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "task/include/task.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RadixSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void RadixSort(std::vector<int> &data);
};

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
