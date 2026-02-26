#pragma once

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "task/include/task.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

class RadixSortSimpleMergingSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RadixSortSimpleMergingSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging